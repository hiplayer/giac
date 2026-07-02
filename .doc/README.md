# Giac 测试功能文档

本目录根据当前测试用例（`bin/` 脚本 + `giac-1.5.0/check/` 回归套件）整理，
描述**已有测试覆盖的 CAS 功能**、对应源码模块与 CTest 名称。

## 文档索引

| 文件 | 内容 |
|------|------|
| [functional-coverage.md](functional-coverage.md) | 按功能域划分的 API 覆盖与源码模块映射 |
| [test-inventory.md](test-inventory.md) | 全部 51 个测试输入文件清单 |
| [gaps.md](gaps.md) | 尚未覆盖或覆盖薄弱的功能域 |
| [rust-migration-plan.md](rust-migration-plan.md) | Giac → Rust 主迁移方案 |
| [cas-long-term-vision.md](cas-long-term-vision.md) | CAS 长远愿景（数学边界、能力分层、参数化语义） |
| [rust-migration-supplement.md](rust-migration-supplement.md) | 移植补充分析（依赖、文法、功能全景） |
| [parser-token-map.md](parser-token-map.md) | Token / 优先级 / 多模式 / Context 映射 |
| [builtin-api-map.md](builtin-api-map.md) | 测试 API → `at_*` → Rust crate |
| [module-division.md](module-division.md) | C++ 算法源文件分工 |
| [algorithm-expr-api.md](algorithm-expr-api.md) | **通用** 稳定/临时 API、Expr↔Poly 边界、防形式漂移（Cursor：`algorithm-expr-api.mdc`） |
| [giac-calculus-api-stability.md](giac-calculus-api-stability.md) | **giac-calculus** crate 稳定/Partial/Pipeline 分层与源码注释约定 |
| [limit-engine-expr-api.md](limit-engine-expr-api.md) | limit_engine 子域 I/O 契约（`.doc/` 专项，无 per-module Cursor 规则） |
| [conformance-testing.md](conformance-testing.md) | Golden / assert_equiv；**`cargo test-timeout` 推荐用法** |
| [GIAC-limit-layered-pipeline](issues/GIAC-limit-layered-pipeline.md) | limit 四层管线落地缺口 |
| [GIAC-algorithm-gaps-open](issues/GIAC-algorithm-gaps-open.md) | **giac-rs 算法未实现总览**（`#[ignore]` 中 4 项 + NotImplemented 索引） |
| [GIAC-core-upstream-gaps](issues/GIAC-core-upstream-gaps.md) | **giac-core crate** vs upstream 主要功能缺口（AlgExtC eval / assume / K 上管线 / 类群 / simplify） |
| [external-test-resources.md](external-test-resources.md) | Maxima / SymPy / Rubi 外部 fixture 与抽取脚本 |
| [known-divergences.md](known-divergences.md) | Rust 与 giac 已知偏离登记 |

## 测试体系概览

```
CTest (53)
├── tommath_smoke_test          libtommath 大整数冒烟
├── giac_bin_*  (41)            bin/ 覆盖驱动脚本（无 golden diff）
├── giac_check_* (10)           check/ 回归测试（stdout diff golden .out）
└── giac_check_factor_segfault_known   已知 segfault（WILL_FAIL）
```

### 两类测试

**覆盖驱动（`bin/test_*`）** — 执行 giac 表达式，驱动 gcov 覆盖，不校验输出正确性。

**回归测试（`check/test*`）** — 将 giac stdout 与黄金文件 `.out` 做 `diff`；
`cas` 和 `geo` 经 `cmake/giac_check_normalize.sh` 归一化浮点/绘图 ID。

## 统计摘要

- 测试输入文件：**51** 个（bin **41** + check **10**）
- CTest 目标：**53** 个（41 bin + 10 check + 1 tommath + 1 segfault）
- 表达式总数：约 **1210** 条
- 涉及 API：**~210** 种 **MVP API**（**M**inimum **V**iable **P**roduct **API**；giac 全库 ~1852 个 `at_*` 注册）— 定义见 [builtin-api-map.md](builtin-api-map.md) §「MVP API 是什么」
