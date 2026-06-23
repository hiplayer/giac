# P3-6 — `poly_algext_roots` 算法规格

**状态:** **resolved**（normative §0–§11 已落地并归档）  
**归档:** [issues_resolved/GIAC-poly-p3-6-roots-algorithm-spec.md](../issues_resolved/GIAC-poly-p3-6-roots-algorithm-spec.md)  
**快照:** 2026-06-23

---

## 实现程度（摘要）

| 范围 | 状态 | 证据 |
|------|------|------|
| deg 1–4 `poly_algext_roots` | ✅ | `roots_dispatch` → linear / quadratic / `split_monic_cubic` / `quartic_roots` |
| 三次统一分裂 C1/C2 | ✅ | `split_monic_cubic_roots_in_session`；纯三次 ω；p≠0 deflate+F1 |
| 四次 F2 resolvent + F3 Euler | ✅ | `euler_four_roots_vanish`、`roots_quartic_t4_plus_t_plus_1` |
| F1 `sqrt_in_field` | ✅ | 塔扫描 + 有界枚举 interim |
| 维数门禁 | ✅ | `D_HARD=24`、`D_RESOLVENT=6`；`t⁴+t+1` baseline dim=24 |
| deg ≥ 5 | N/A | `NotImplemented`（solve S0 单支 `rootof`） |

```bash
cargo test -p giac-core poly_roots::tests --release   # 25 passed, 0 ignored
```

**缺口索引（仍 open）：** [GIAC-poly-p3-6-quartic-roots-gaps.md](GIAC-poly-p3-6-quartic-roots-gaps.md)（F4′ / F5 plan、solve S7）

---

## §12 升级路径（plan，仍 open）

**不阻塞 P3-6 验收。** 完整 normative 背景见 [归档 §12](../issues_resolved/GIAC-poly-p3-6-roots-algorithm-spec.md)。

**现状（interim）：** F1 `try_sqrt_in_field` 用塔层扫描 + 有界 ±1 枚举；Euler 第二开方用 shallow + blind adjoin。根正确（`verify_root`），但 S₄ 四次（`t⁴+t+1`）实测 `d_L=24`，枚举非长期方案。

| 优先级 | 项 | 替代对象 | 落点 |
|--------|-----|----------|------|
| **P1** | **F5 — 结构开方** | §6 的 1c–1e ±1 枚举 | `ext_tower::try_square_root_in_field` |
| **P2** | **F4′ — Galois 第二 Euler 开方** | `sqrt_in_field_euler_second` shallow+adjoin | `galois_automorphism.rs` + `poly_roots.rs` |
| P3 | F5 收敛 | `FieldSession::try_sqrt_in_field` 薄封装 | `field_session.rs` |

**依赖顺序：** F1 ✅ → **F5** → **F4′** → 删除 `TRY_SQRT_*` interim 常量。

**plan DoD：**

- [ ] F4′：Gal≅A₄ **已证** 四次上 \(d_L \le D_QUARTIC_A4_REF\) 且四根 verify；或登记为何 \(d_L>12\)
- [ ] `t⁴+t+1`：保持 \(d_L \le D_HARD\) + 四根 verify（现行 baseline）
- [ ] `try_square_root` 无 dim³ 枚举热路径；`poly_roots::tests` 仍 ≤ `T_ROOTS`
- [ ] [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md) §F5 / F4′ 同步

**跟踪 issue：** [GIAC-poly-quartic-roots-F1-F5](GIAC-poly-quartic-roots-F1-F5.md)（F5 🟡、F4′ open）
