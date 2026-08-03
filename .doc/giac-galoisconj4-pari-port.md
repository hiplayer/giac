# GIAC `galoisconj4` — Pari `galconj.c` 函数对照

**上游：** `pari/src/basemath/galconj.c`（基线 giac-rs 对标 Pari，非 giac C++ 行级）  
**Rust：** `giac-core/src/algebra/galoisconj4/` + `algebra/galois_conj.rs`  
**测试规格：** [GIAC-galoisconj4-port-plan.md](issues/GIAC-galoisconj4-port-plan.md) G1–G6；shadow 见 `upstream_shadow` 模块

**说明：** `galoisconj4` = p-adic Frobenius 提升 + 模 Vandermonde + `permtopol`，**不是** `buch2.c` LLL。

---

## Field 层入口（`galois_conj.rs`）

| Pari | Rust | 路由 |
|------|------|------|
| `galoisconj_monic(T)` | `galoisconj_upstream_only` / `FieldGaloisSnapshot::compute` | deg 1/2 → `[α]` / `[−α,α]`；`numberofconjugates==1` → `[α]`；else `galois_init` → g1 |
| `galoisconj(nf)` | `galoisconj_in_field` → `GaloisConjugates::compute` | 同上（经 snapshot） |
| `galoisinit(T, flag)` | `galois_init(t, flag)` | 一次 g4；P6 缓存在 `FieldGaloisSnapshot` |
| `nfgaloismatrix(nf, σ)` | `nfgaloismatrix` / `NfAutMatrix` | 幂基 ℤ-矩阵 |
| `automorphism_matrices` | `AutomorphismMatrices::compute` | 读 snapshot |
| `FB_aut_perm` | `FbAutPerms` | 因子基 ideal 置换 |
| `automorphism_perms` | `EmbAutPerms::from_galois_init` | g4 成功；失败 → arch 启发式（BNF 专用） |
| grow `idealperm`+`embperm` | `GaloisAutPerms` / `GrhRelCache.galois_snapshot` | P6 同源 |

**已移除（2026-08-03 P5-upstream）：** Rust arch easy 主路径（`galoisconj_easy`、f64 `sigma_alpha_from_arch_slot_perm`）。`galois_root_perms` 仅服务 `EmbAutPerms` 回退，**不参与** conjugates。

---

## `galoisconj4_main` 管线（`galoisconj4/`）

| Pari | Rust 模块 | 备注 |
|------|-----------|------|
| `galoisanalysis` | `analysis.rs` | WSS / 可表群判定 |
| `numberofconjugates` | `analysis::expected_conjugate_count` | G1 计数 |
| `galoisborne` | `borne.rs` | p-adic 精度界 |
| `galois_find_frobenius` | `frobenius.rs` | 素数扫描 + Frobenius |
| `galoisgen` / `galoisgenlift` | `gen.rs`, `lift.rs` | 生成元提升 |
| `galoisgenfixedfield` | `fixed_field.rs` | 固定域子问题 |
| A₄ / S₄ / F₃₆ 快路 | `specials.rs` | P3c |
| `vectopol` / `permtopol` | `perm.rs` | P4；`galois_vec_perm_to_pol` |
| `galoisconj4_main` | `main.rs::galoisconj4_main` | 编排 |
| `galoisconj1` | `galoisconj1.rs` | P7 回退 |
| `nfroots` | `nf_roots`（经 g1） | P7 |

---

## 测试与 CI

| 层级 | 谓词 | Rust 锚点 |
|------|------|-----------|
| G1–G4 | 计数 / 最小多项式 / 自同构 / identity-last | `galoisconj_probe_*_g1_g4` |
| G6 | Pari nfelt multiset | `g6_pari_nfelt::*` |
| Shadow | entry ≡ upstream-only ≡ snapshot | `upstream_shadow::*` |
| Golden 脚本 | 探针批跑 | `giac-rs/scripts/galoisconj_golden.sh` |

```bash
cd giac-rs
./scripts/galoisconj_golden.sh
PARI_GOLDEN=1 ./scripts/galoisconj_golden.sh   # 可选 gp 对照
```

---

## 已知余量

- **P7：** 非 WSS、`1<c<n` 探针 G6 金值待补
- **P8：** `pr_orbit_fill` → `be_honest`
- **EmbAutPerms：** g1-only（无 `galoisinit`）仍 arch 启发式 — 见 [known-divergences.md](known-divergences.md) DIV-105
