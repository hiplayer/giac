//! GIAC-202: `bin/test_desolve*` parse smoke.

use giac_conformance::upstream_root;
use giac_ode::xcas_default;
use giac_parse::parse_program;

#[test]
fn test_desolve_parse_full_file() -> Result<(), String> {
    let path = upstream_root().join("bin/test_desolve");
    let input = std::fs::read_to_string(&path)
        .map_err(|e| format!("read {}: {e}", path.display()))?;
    parse_program(&input, &xcas_default()).map_err(|e| format!("parse test_desolve: {e}"))?;
    Ok(())
}

#[test]
fn test_desolve_ext_parse_full_file() -> Result<(), String> {
    let path = upstream_root().join("bin/test_desolve_ext");
    let input = std::fs::read_to_string(&path)
        .map_err(|e| format!("read {}: {e}", path.display()))?;
    parse_program(&input, &xcas_default()).map_err(|e| format!("parse test_desolve_ext: {e}"))?;
    Ok(())
}
