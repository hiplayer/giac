use std::io::Write;
use std::process::{Command, Stdio};

fn run_stdin(input: &str) -> String {
    let mut child = Command::new(env!("CARGO_BIN_EXE_giac-cli"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn giac-cli");
    write!(child.stdin.as_mut().expect("stdin"), "{input}").expect("write stdin");
    drop(child.stdin.take());
    let output = child.wait_with_output().expect("wait giac-cli");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout).trim().to_string()
}

#[test]
fn cli_evaluates_from_stdin() {
    assert_eq!(run_stdin("gcd(45,75);"), "15");
}

#[test]
fn cli_multiple_statements() {
    assert_eq!(run_stdin("1+1;2+2;"), "2,4");
}

#[test]
fn cli_assignment_prints_value() {
    assert_eq!(run_stdin("x:=5;"), "5");
}

#[test]
fn cli_reads_script_file() {
    let path = std::env::temp_dir().join("giac_cli_test_script.giac");
    std::fs::write(&path, "normal((x+1)^2);").expect("write temp script");
    let output = Command::new(env!("CARGO_BIN_EXE_giac-cli"))
        .arg(&path)
        .output()
        .expect("run giac-cli with file");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(stdout.trim(), "x^2+2*x+1");
    let _ = std::fs::remove_file(path);
}
