use unicode_normalization::UnicodeNormalization;

pub fn normalize_barcode(raw: &str, idx: usize) -> String {
    let cleaned = raw.trim().nfc().collect::<String>();
    if cleaned.is_empty() {
        return format!("cell_{:08}", idx + 1);
    }
    cleaned
}

pub fn normalize_gene_symbol(gene_id: &str, gene_symbol: Option<&str>, idx: usize) -> String {
    let candidate = gene_symbol.unwrap_or(gene_id).trim();
    let mut value = if candidate.is_empty() {
        format!("gene_{:08}", idx + 1)
    } else {
        candidate.nfc().collect::<String>()
    };

    if (value.starts_with("ENSG") || value.starts_with("ENSMUSG"))
        && value
            .rsplit_once('.')
            .map(|(_, s)| s.chars().all(|c| c.is_ascii_digit()))
            .unwrap_or(false)
    {
        value = value.split('.').next().unwrap_or(&value).to_string();
    }

    value
}
