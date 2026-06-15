//! GIAC-203: `bin/test_integrate*` parse smoke.

use giac_conformance::upstream_root;
use giac_ode::xcas_default;
use giac_parse::parse_program;

#[test]
fn test_integrate_parse_full_file() -> Result<(), String> {
    let path = upstream_root().join("bin/test_integrate");
    let input = std::fs::read_to_string(&path)
        .map_err(|e| format!("read {}: {e}", path.display()))?;
    parse_program(&input, &xcas_default()).map_err(|e| format!("parse test_integrate: {e}"))?;
    Ok(())
}

#[test]
fn test_integrate_more_parse_full_file() -> Result<(), String> {
    let path = upstream_root().join("bin/test_integrate_more");
    let input = std::fs::read_to_string(&path)
        .map_err(|e| format!("read {}: {e}", path.display()))?;
    parse_program(&input, &xcas_default())
        .map_err(|e| format!("parse test_integrate_more: {e}"))?;
    Ok(())
}

#[test]
fn test_integrate_ext_parse_full_file() -> Result<(), String> {
    let path = upstream_root().join("bin/test_integrate_ext");
    let input = std::fs::read_to_string(&path)
        .map_err(|e| format!("read {}: {e}", path.display()))?;
    parse_program(&input, &xcas_default()).map_err(|e| format!("parse test_integrate_ext: {e}"))?;
    Ok(())
}
