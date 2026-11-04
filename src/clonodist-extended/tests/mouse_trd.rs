use clonodist_extended::chain_config;

#[test]
fn mouse_delta_requires_explicit_reference() {
    assert!(chain_config("delta", "mouse").is_err());
    assert!(chain_config("gamma,delta", "mouse").is_err());
    assert!(chain_config("gamma", "mouse").is_ok());
}
