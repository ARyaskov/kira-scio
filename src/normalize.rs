//! Canonicalization helpers for gene/barcode strings.

use unicode_normalization::UnicodeNormalization;

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
    let cleaned = if candidate.is_empty() {
        synth_gene(idx)
    } else {
        nfc(candidate)
    };
    strip_ensembl_version(&cleaned)
}

/// Prefers `gene_id` over `fallback_symbol`, then synthesises a placeholder.
pub fn normalize_gene_id(gene_id: &str, fallback_symbol: Option<&str>, idx: usize) -> String {
    let id_trimmed = gene_id.trim();
    let value = if !id_trimmed.is_empty() {
        nfc(id_trimmed)
    } else {
        match fallback_symbol.map(str::trim).filter(|s| !s.is_empty()) {
            Some(s) => nfc(s),
            None => synth_gene(idx),
        }
    };
    strip_ensembl_version(&value)
}

fn nfc(value: &str) -> String {
    value.nfc().collect()
}

fn nfc_trim(value: &str) -> String {
    value.trim().nfc().collect()
}

fn strip_ensembl_version(value: &str) -> String {
    if (value.starts_with("ENSG") || value.starts_with("ENSMUSG"))
        && value
            .rsplit_once('.')
            .map(|(_, s)| !s.is_empty() && s.chars().all(|c| c.is_ascii_digit()))
            .unwrap_or(false)
    {
        return value.split('.').next().unwrap_or(value).to_string();
    }
    value.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn ensembl_version_is_stripped() {
        assert_eq!(normalize_gene_id("ENSG00000001.5", None, 0), "ENSG00000001");
        assert_eq!(
            normalize_gene_id("ENSMUSG00000000001.2", None, 0),
            "ENSMUSG00000000001"
        );
        // Non-versioned ids should pass through.
        assert_eq!(normalize_gene_id("ENSG00000001", None, 0), "ENSG00000001");
    }
}
