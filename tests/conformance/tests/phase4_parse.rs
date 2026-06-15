//! Phase 4 bin scripts: parse-only smoke (GIAC-201–204).

use giac_conformance::{script_lines, upstream_root, PHASE4_SCRIPTS};
use giac_ode::xcas_default;
use giac_parse::parse_program;

#[test]
fn phase4_scripts_parse() -> Result<(), String> {
    let ctx = xcas_default();
    for name in PHASE4_SCRIPTS {
        let path = upstream_root().join("bin").join(name);
        let input = std::fs::read_to_string(&path)
            .map_err(|e| format!("read {}: {e}", path.display()))?;
        parse_program(&input, &ctx).map_err(|e| format!("parse {name}: {e}"))?;
        let lines = script_lines(&path)?;
        assert!(!lines.is_empty(), "{name} has no lines");
    }
    Ok(())
}
