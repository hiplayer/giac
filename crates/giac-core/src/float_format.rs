/// Format a float with `digits` significant figures (giac `evalf` style).
pub fn format_float(v: f64, digits: u32) -> String {
    if v.is_nan() {
        return "undef".to_string();
    }
    if v.is_infinite() {
        return if v.is_sign_positive() {
            "inf".to_string()
        } else {
            "-inf".to_string()
        };
    }
    if v == 0.0 {
        return "0".to_string();
    }
    let exp = v.abs().log10().floor() as i32;
    let scale = 10_f64.powi(digits as i32 - 1 - exp);
    let rounded = (v * scale).round() / scale;
    let s = format!("{rounded:.12}");
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}
