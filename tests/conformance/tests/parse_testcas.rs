use giac_core::Context;
use giac_parse::parse_program;

#[test]
fn parse_testcas_up_to_50() {
    let text = std::fs::read_to_string(
        "/home/kanli.hu/upstream/giac/giac/giac-1.5.0/check/testcas",
    )
    .expect("testcas");
    let lines: Vec<&str> = text.lines().collect();
    let ctx = Context::xcas_default();
    for n in 1..=50 {
        let input: String = lines
            .iter()
            .take(n)
            .map(|l| if l.ends_with(';') { l.to_string() } else { format!("{l};") })
            .collect();
        if let Err(e) = parse_program(&input, &ctx) {
            panic!("parse failed at line {n} ({:?}): {e}", lines.get(n - 1));
        }
    }
}
