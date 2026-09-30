//! Canonicalization helpers for gene/barcode strings.

use unicode_normalization::UnicodeNormalization;

/// UTF-8 byte order mark as it appears at the start of a decoded line.
pub const UTF8_BOM: char = '\u{feff}';

/// Strips a leading UTF-8 BOM (written by Excel and some Windows tools).
/// Apply to the first line of a text file only.
pub fn strip_bom(line: &str) -> &str {
    line.strip_prefix(UTF8_BOM).unwrap_or(line)
}

pub fn normalize_barcode(raw: &str, idx: usize) -> String {
    let cleaned = nfc_trim(raw);
    if cleaned.is_empty() {
        return synth_barcode(idx);
    }
    cleaned
}

pub fn synth_barcode(idx: usize) -> String {
    synth_padded("cell_", idx + 1)
}

pub fn synth_gene(idx: usize) -> String {
    synth_padded("gene_", idx + 1)
}

fn synth_padded(prefix: &'static str, value: usize) -> String {
    const WIDTH: usize = 8;
    let mut buf = itoa::Buffer::new();
    let s = buf.format(value).as_bytes();
    let n = s.len();
    if n >= WIDTH {
        let mut out = String::with_capacity(prefix.len() + n);
        out.push_str(prefix);
        out.push_str(std::str::from_utf8(s).unwrap());
        return out;
    }
    let pad = WIDTH - n;
    let mut out = String::with_capacity(prefix.len() + WIDTH);
    out.push_str(prefix);
    for _ in 0..pad {
        out.push('0');
    }
    out.push_str(std::str::from_utf8(s).unwrap());
    out
}

/// Prefers `gene_symbol` over `gene_id`. Used when only one display name is
/// needed.
pub fn normalize_gene_symbol(gene_id: &str, gene_symbol: Option<&str>, idx: usize) -> String {
    let candidate = gene_symbol.unwrap_or(gene_id).trim();
    if candidate.is_empty() {
        synth_gene(idx)
    } else {
        nfc(candidate)
    }
}

/// Prefers `gene_id` over `fallback_symbol`, then synthesises a placeholder.
pub fn normalize_gene_id(gene_id: &str, fallback_symbol: Option<&str>, idx: usize) -> String {
    let id_trimmed = gene_id.trim();
    if !id_trimmed.is_empty() {
        nfc(id_trimmed)
    } else {
        match fallback_symbol.map(str::trim).filter(|s| !s.is_empty()) {
            Some(s) => nfc(s),
            None => synth_gene(idx),
        }
    }
}

fn nfc(value: &str) -> String {
    value.nfc().collect()
}

fn nfc_trim(value: &str) -> String {
    value.trim().nfc().collect()
}

/// Strips the `.N` version suffix from an Ensembl stable id of any species
/// and feature type (`ENSG…`, `ENSMUSG…`, `ENSDARG…`, `ENST…`, `ENSP…`, …).
/// Ensembl stable ids are `ENS`, an optional species code of up to five
/// upper-case letters, a feature-type code, and at least six digits;
/// anything else is returned unchanged. Applied by the `Reader` unless
/// `ReaderOptions::strip_ensembl_versions` is off.
pub fn strip_ensembl_version(value: &str) -> &str {
    let Some((stem, version)) = value.rsplit_once('.') else {
        return value;
    };
    if version.is_empty() || !version.bytes().all(|b| b.is_ascii_digit()) {
        return value;
    }
    let Some(rest) = stem.strip_prefix("ENS") else {
        return value;
    };
    let letters = rest.bytes().take_while(|b| b.is_ascii_uppercase()).count();
    let digits = rest[letters..].len();
    let all_digits = rest[letters..].bytes().all(|b| b.is_ascii_digit());
    // Species code (0-5 letters) plus at least one feature-type letter.
    if (1..=6).contains(&letters) && digits >= 6 && all_digits {
        stem
    } else {
        value
    }
}

/// Makes labels unique in place the way scanpy's `var_names_make_unique`
/// does: the second occurrence of `X` becomes `X-1`, the third `X-2`, and so
/// on. Returns the number of labels renamed.
pub fn make_labels_unique(labels: &mut [String]) -> usize {
    use rustc_hash::FxHashMap;
    let mut seen: FxHashMap<String, usize> = FxHashMap::default();
    for l in labels.iter() {
        *seen.entry(l.clone()).or_insert(0) += 1;
    }
    let mut counters: FxHashMap<String, usize> = FxHashMap::default();
    let mut renamed = 0;
    for l in labels.iter_mut() {
        if seen.get(l).copied().unwrap_or(0) <= 1 {
            continue;
        }
        let n = counters.entry(l.clone()).or_insert(0);
        if *n > 0 {
            let mut candidate = format!("{l}-{n}");
            // Avoid colliding with a label that already exists in the input.
            while seen.contains_key(&candidate) {
                *n += 1;
                candidate = format!("{l}-{n}");
            }
            *l = candidate;
            renamed += 1;
        }
        *n += 1;
    }
    renamed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strip_bom_removes_only_a_leading_mark() {
        assert_eq!(strip_bom("\u{feff}gene\tC1"), "gene\tC1");
        assert_eq!(strip_bom("gene\tC1"), "gene\tC1");
        assert_eq!(strip_bom("a\u{feff}b"), "a\u{feff}b");
        assert_eq!(strip_bom(""), "");
    }

    #[test]
    fn synth_barcode_is_zero_padded() {
        assert_eq!(synth_barcode(0), "cell_00000001");
        assert_eq!(synth_barcode(9), "cell_00000010");
        assert_eq!(synth_barcode(99_999_999), "cell_100000000");
    }

    #[test]
    fn synth_gene_format_matches_legacy() {
        assert_eq!(synth_gene(0), "gene_00000001");
        assert_eq!(synth_gene(7), "gene_00000008");
    }

    #[test]
    fn normalize_id_prefers_id_then_symbol() {
        assert_eq!(normalize_gene_id("ENSG1", Some("BRCA1"), 0), "ENSG1");
        assert_eq!(normalize_gene_id("", Some("BRCA1"), 0), "BRCA1");
        assert_eq!(normalize_gene_id("", None, 4), "gene_00000005");
    }

    #[test]
    fn normalize_symbol_prefers_symbol_then_id() {
        assert_eq!(normalize_gene_symbol("ENSG1", Some("BRCA1"), 0), "BRCA1");
        assert_eq!(normalize_gene_symbol("ENSG1", None, 0), "ENSG1");
        assert_eq!(normalize_gene_symbol("", None, 4), "gene_00000005");
    }

    #[test]
    fn ensembl_version_is_stripped_for_any_species_and_feature_type() {
        for (input, expected) in [
            ("ENSG00000139618.15", "ENSG00000139618"),
            ("ENSMUSG00000000001.2", "ENSMUSG00000000001"),
            ("ENSDARG00000000001.3", "ENSDARG00000000001"),
            ("ENSRNOG00000000002.1", "ENSRNOG00000000002"),
            ("ENSSSCG00000000003.7", "ENSSSCG00000000003"),
            ("ENST00000380152.8", "ENST00000380152"),
            ("ENSP00000369497.3", "ENSP00000369497"),
            // Unversioned or non-Ensembl labels pass through.
            ("ENSG00000139618", "ENSG00000139618"),
            ("BRCA2", "BRCA2"),
            ("MT-CO1", "MT-CO1"),
            ("FBgn0000001.1", "FBgn0000001.1"),
            ("ENS.5", "ENS.5"),
            ("ENSG1.5", "ENSG1.5"),
            ("ENSG00000139618.x", "ENSG00000139618.x"),
        ] {
            assert_eq!(strip_ensembl_version(input), expected, "{input}");
        }
        // Normalization itself no longer strips; the Reader applies it.
        assert_eq!(
            normalize_gene_id("ENSG00000001.5", None, 0),
            "ENSG00000001.5"
        );
    }

    #[test]
    fn make_labels_unique_appends_running_suffixes() {
        let mut v: Vec<String> = ["A", "B", "A", "A", "B-1", "B"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let renamed = make_labels_unique(&mut v);
        assert_eq!(renamed, 3);
        assert_eq!(v, vec!["A", "B", "A-1", "A-2", "B-1", "B-2"]);
        let mut unique: Vec<String> = vec!["x".into(), "y".into()];
        assert_eq!(make_labels_unique(&mut unique), 0);
    }
}
