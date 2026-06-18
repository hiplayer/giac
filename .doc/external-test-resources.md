# 外部开源测试资源

**状态:** normative（fixture 与抽取脚本规格）  
**相关:** [conformance-testing.md](conformance-testing.md)、[phase4-issues.md](phase4-issues.md) §2.6、[cas-long-term-vision.md](cas-long-term-vision.md) Phase B

本文定义 **Maxima / SymPy / Rubi** 等外部语料如何进入 giac-rs conformance：**目录布局、JSON schema、抽取脚本、验证层级、CI 门禁**。

giac 自有 golden（53 CTest）仍为 **T5 交付门禁**；外部资源为 **T1 扩展规格**，不替代 `check/*.out` 字面 diff。

---

## 1. 验证层级（与 phase4 对齐）

| 层级 | 工具 | 角色 | PR 必过？ |
|------|------|------|-----------|
| **T5** | giac `bin/` + `check/` | 回归 golden | ✅ MVP |
| **T1** | `sympy_verify.py` | 数学 oracle | ✅ 启用 fixture 行 |
| **T2** | Rust `assert_equiv` | 等价不同形 | 可选回退 |
| **T3** | `proptest` | 属性（gcd、diff∘integrate） | 可选 |
| **T4** | giac C++ 双跑 | 对照参考实现 | 可选 |
| **T-ext** | 外部 fixture JSON | 扩展覆盖 | 仅 `enabled: true` 行 |

**合并规则：** 外部 fixture 行启用后，须通过 **T1**（或该条 `verify` 字段声明的模式）；不比 Maxima/SymPy **字面**输出。

---

## 2. 目录结构（giac-rs）

根路径：`giac-rs/tests/conformance/`（与 [rust-migration-plan.md](rust-migration-plan.md) conformance harness 并列）。

```
giac-rs/tests/conformance/
├── fixtures/
│   ├── README.md                 # 本布局说明（副本见下文 §3）
│   ├── _schema/
│   │   └── fixture-v1.schema.json
│   ├── internal/                 # 文档索引，不复制 giac golden
│   │   └── README.md
│   ├── maxima/                   # Maxima rtest 抽取结果
│   │   ├── limit.draft.json
│   │   ├── limit.json            # 人工审核后（去 .draft）
│   │   ├── integrate.draft.json
│   │   ├── taylor.draft.json
│   │   └── gruntz.draft.json     # share/tests/rtest_gruntz.mac
│   ├── sympy/                    # SymPy 测试抽取 / 手工 curated
│   │   ├── limits.draft.json
│   │   ├── integrals.draft.json
│   │   └── polys.draft.json
│   └── rubi/                     # Rubi 子集（抽样，非全量）
│       ├── manifests/
│       │   └── catalog.yaml      # 源文件 → fixture 映射
│       └── samples/
│           ├── timofeev_705.draft.json
│           └── trig_4_3_11.draft.json
├── sources/                      # 上游克隆（.gitignore，本地可选）
│   ├── README.md
│   ├── maxima/                   # git clone maxima
│   ├── sympy/
│   └── rubi/
├── scripts/
│   ├── extract_maxima_rtest.py   # Maxima .mac → JSON（§5）
│   ├── extract_sympy_tests.py    # SymPy test_*.py → JSON（§6，待实现）
│   ├── extract_rubi_problems.py  # Rubi 纯文本 → JSON（§7，待实现）
│   └── sympy_verify.py           # T1 oracle（见 phase4-issues §2.6.2）
└── tests/
    ├── phase4_maxima_rtest.rs    # 兼容旧名；读 fixtures/maxima/limit.json
    ├── external_maxima.rs        # 按域聚合 maxima/*.json
    ├── external_sympy.rs
    └── external_rubi.rs          # nightly 或 `#[ignore]` 默认
```

### 2.1 命名约定

| 后缀 | 含义 |
|------|------|
| `*.draft.json` | 脚本生成，**默认全部 `enabled: false`**，待人工审核 |
| `*.json` | 审核通过，可部分 `enabled: true` 进 CI |
| `id` 前缀 | `LIM-M-` Maxima limit · `INT-M-` Maxima integrate · `LIM-S-` SymPy · `INT-R-` Rubi |

### 2.2 与旧路径兼容

| 旧路径（phase4） | 新路径 |
|------------------|--------|
| `fixtures/phase4_maxima_rtest.json` | `fixtures/maxima/limit.json`（迁移中；旧格式见 §2.3） |
| `fixtures/phase4_integrate_table.json` | 保留；属 **内部 curated**，非外部抽取 |

### 2.3 旧格式兼容（`phase4_maxima_rtest.json`）

现有 `phase4_maxima_rtest.rs` 读取的 JSON 为 **legacy 子集**，每条 entry **必须**含 `api` 字段：

```json
{
  "version": 1,
  "entries": [
    {
      "id": "MX_LIMIT_WESTER-LIM-001",
      "line": "limit((1 + 1/n)^n, n, +infinity)",
      "api": "limit",
      "verify": "sympy_limit",
      "enabled": true,
      "maxima_input": "...",
      "maxima_expected": "..."
    }
  ]
}
```

`extract_maxima_rtest.py` v1 输出 **fixture-v1** 全字段，并自动填充 legacy 别名（`api`、`maxima_input`、`maxima_expected`、`source_file`、`source_line`）。迁移路径：

1. 新抽取 → `fixtures/maxima/*.draft.json`
2. 审核后合并或替换 `phase4_maxima_rtest.json`
3. 长期：`external_maxima.rs` 聚合 `fixtures/maxima/*.json`

---

## 3. JSON Fixture 规格（`fixture-v1`）

Schema：`fixtures/_schema/fixture-v1.schema.json`。

### 3.1 顶层

```json
{
  "version": 1,
  "source": {
    "name": "maxima",
    "file": "tests/rtest_limit.mac",
    "revision": "optional-git-rev",
    "url": "https://sourceforge.net/p/maxima/code/ci/master/tree/tests/rtest_limit.mac",
    "license": "GPL-2.0-or-later",
    "extracted_at": "2026-06-18",
    "extractor": "extract_maxima_rtest.py",
    "extractor_version": "1"
  },
  "domain": "limit",
  "verify": "sympy_limit",
  "entries": []
}
```

### 3.2 `domain` 枚举

`limit` · `integrate` · `diff` · `solve` · `taylor` · `series` · `partfrac` · `desolve` · `poly` · `trig` · `other`

### 3.3 `verify` 枚举（对应 `sympy_verify.py` 子命令）

| `verify` | 判定 |
|----------|------|
| `sympy_limit` | `sp.limit(f, x, a)` |
| `integrate_derivative` | `simplify(diff(F,x)-f)==0` |
| `diff_inverse` | `simplify(diff(f,x)-F)==0` |
| `solve_residual` | 每解代入方程残差为 0 |
| `partfrac_expand` | 部分分式展开等于原式 |
| `desolve_odesol` | `checkodesol` |
| `parse_only` | 仅 `giac-parse`，不 eval |
| `literal` | 与 `expected` 字面比（少用） |

### 3.4 `entries[]` 字段

| 字段 | 必填 | 说明 |
|------|------|------|
| `id` | ✅ | 全局唯一，如 `LIM-M-014` |
| `line` | ✅ | **giac 语法**单行：`limit(sin(x)/x,x,0)` |
| `enabled` | ✅ | 默认 `false`（draft） |
| `source_ref` | 推荐 | 上游定位：`rtest_limit.mac:42` |
| `upstream_input` | 可选 | 抽取前原文（Maxima/SymPy） |
| `expected_hint` | 可选 | 人读参考；**不作 CI 断言** |
| `tags` | 可选 | `gruntz`, `+infinity`, `CK-INT-58`, … |
| `giac_issue` | 可选 | 启用门槛，如 `GIAC-215` |
| `bin_ref` | 可选 | 关联 giac `bin/test_*` 行 |
| `notes` | 可选 | 审核备注 |

---

## 4. 上游源清单

### 4.1 Maxima（GPL-2.0+）

| 文件 | 域 | 优先级 |
|------|-----|--------|
| `tests/rtest_limit.mac`, `rtest_limit_extra.mac` | limit | P0（已部分抽取） |
| `share/` 下 `rtest_gruntz.mac` | limit / gruntz | P0 |
| `tests/rtest_integrate.mac`, `rtestint.mac` | integrate | P1 |
| `tests/rtest_taylor.mac`, `rtest_powerseries.mac` | series | P1 |
| `tests/rtest_trig.mac` | trig | P2 |
| `tests/rtest_gcd.mac` | poly | P2 |
| `tests/rtestode.mac` | desolve | P2 |
| `tests/wester_problems/` | 综合 | P3 |

克隆：`sources/maxima` → [Maxima git](https://sourceforge.net/p/maxima/code/ci/master/tree/)

### 4.2 SymPy（BSD-3-Clause）

按模块抽取 `sympy/*/tests/test_*.py` 中带明确断言的用例：

| 路径 | 域 |
|------|-----|
| `sympy/series/tests/test_limits.py` | limit |
| `sympy/integrals/tests/test_integrals.py` | integrate |
| `sympy/polys/tests/test_gcdex.py` 等 | poly |
| `sympy/solvers/tests/test_solvers.py` | solve |

克隆：`sources/sympy` → https://github.com/sympy/sympy

**策略：** 首版 **手工 curated** + `extract_sympy_tests.py` 辅助；输出 `fixtures/sympy/*.draft.json`。

### 4.3 Rubi（MIT）

| 资源 | 规模 | 用法 |
|------|------|------|
| [rulebasedintegration.org/testProblems.html](https://rulebasedintegration.org/testProblems.html) | ~72k | 按目录抽样 |
| Timofeev 独立集 | 705 | `rubi/samples/timofeev_705.draft.json` |
| Nasser Abbasi 2024 报告 | 交叉基准 | 人工对照，不自动导入 |

格式：`(integrand, variable)` 纯文本 → giac `integrate(...,x)`；验证 **`integrate_derivative`**。

**CI：** 默认 `enabled: false` 或单独 `cargo test --test external_rubi -- --ignored`；全量仅 nightly。

### 4.4 其他（文档索引）

| 资源 | 许可 | 见 |
|------|------|-----|
| FriCAS `src/input/*.input` | BSD | `sources/README.md` |
| Wester benchmark | 公有题干 | Maxima `wester_problems` |
| giac `bin/` + `check/` | GPL | [conformance-testing.md](conformance-testing.md) |

---

## 5. `extract_maxima_rtest.py` 规范

**路径：** `giac-rs/tests/conformance/scripts/extract_maxima_rtest.py`

### 5.1 职责

1. 读取 Maxima `tests/*.mac`（及 `share/**/rtest_gruntz.mac`）
2. 识别 `limit` / `integrate` / `diff` / `solve` / `taylor` 调用与期望结果
3. Maxima → giac 语法转换（§5.3）
4. 输出 **`*.draft.json`**，`enabled: false`

### 5.2 CLI

```bash
cd giac-rs/tests/conformance

# 单文件
python3 scripts/extract_maxima_rtest.py \
  --input sources/maxima/tests/rtest_limit.mac \
  --domain limit \
  --output fixtures/maxima/limit.draft.json \
  --prefix LIM-M

# 批量（manifest）
python3 scripts/extract_maxima_rtest.py \
  --manifest scripts/manifests/maxima.yaml

# 预览，不写文件
python3 scripts/extract_maxima_rtest.py \
  --input sources/maxima/tests/rtest_limit.mac \
  --domain limit \
  --dry-run

# 仅抽取含 tag 的用例
python3 scripts/extract_maxima_rtest.py \
  --input sources/maxima/share/gruntz/rtest_gruntz.mac \
  --domain limit --tag gruntz \
  --output fixtures/maxima/gruntz.draft.json
```

| 参数 | 说明 |
|------|------|
| `--input PATH` | 单个 `.mac` 文件（可重复） |
| `--manifest PATH` | YAML：`[{input, domain, output, prefix}]` |
| `--domain NAME` | 覆盖/默认 domain |
| `--output PATH` | 输出 JSON（必填，除非 `--dry-run`） |
| `--prefix STR` | entry `id` 前缀，默认 `MAX-M` |
| `--verify NAME` | 默认 verify；limit→`sympy_limit`，integrate→`integrate_derivative` |
| `--tag STR` | 写入每条 `tags[]` |
| `--start-id N` | 数字后缀起始，默认 1 |
| `--include-known-fail` | 包含 Maxima testsuite 已知失败标记块 |
| `--dry-run` | 打印 JSON 到 stdout |
| `--verbose` | 解析警告 |

退出码：`0` 成功；`1` 参数错误；`2` 未抽取到任何条目。

### 5.3 Maxima → giac 转换（normative）

| Maxima | giac |
|--------|------|
| `%pi`, `%e` | `pi`, `exp(1)` |
| `inf`, `minf`, `infinity` | `+infinity` |
| `-inf`, `-infinity` | `-infinity` |
| `^` | `^`（不变） |
| `*` 显式乘 | 保持或规范化（解析器兼容即可） |
| `limit(f,x,a)` | `limit(f,x,a)` |
| `limit(f,x,inf)` | `limit(f,x,+infinity)` |
| `integrate(f,x)` | `integrate(f,x)` |
| `integrate(f,x,a,b)` | `integrate(f,x,a,b)` |
| `diff(f,x)` | `diff(f,x)` |
| `taylor(f,x,a,n)` | `series(f,x,a,n)` 或 `taylor(...)`（按 domain 配置） |
| `sqrt`, `sin`, `cos`, `log` | 同名（giac 用 `ln` 时脚本将 `log`→`ln`） |
| `true`, `false` | 跳过（非 CAS 表达式） |
| `und`, `ind` | 保留在 `expected_hint`；`line` 仍输出，`tags: ["maxima_special"]` |

**不做：** 完整 Maxima 语法翻译（`block`, `matchdeclare`, `ask` 交互）；此类块 **跳过** 并 `--verbose` 警告。

### 5.4 解析模式（实现要求）

按优先级匹配：

1. **rtest 块：** `/* comment */` 后连续两行 `expr;` / `expected;`
2. **单行 assert：** `[FXX]: expr;` 或 `expr;` 下一行 `expected;`
3. **limit/integrate 行：** 行内 `limit(...)` / `integrate(...)`，期望来自同行 `=` 或下一行
4. **gruntz 命名：** 注释含 `CK-INT-` → `tags` 追加 `CK-INT-xx`

抽取失败行写入 stderr，不中断整文件。

### 5.5 输出与审核流程

```
extract → *.draft.json (enabled: false)
    → 人工 diff line / tags / giac_issue
    → mv limit.draft.json → limit.json，部分 enabled: true
    → cargo test -p giac-conformance --test external_maxima
```

**禁止：** 脚本直接写入 `enabled: true`（除非未来 `--force-enable` 且 PR 明确标注）。

### 5.6 版本与扩展

- `extractor_version` 写入 JSON `source.extractor_version`
- **破坏性变更**（schema、id 规则、转换语义）→ 版本 +1，文档更新 §3
- 新增 domain：同步扩展 `fixture-v1.schema.json`、`sympy_verify.py`、`extract_maxima_rtest.py` 的 `DOMAIN_VERIFY`

---

## 6. `extract_sympy_tests.py`（规格摘要，待实现）

| 项 | 约定 |
|----|------|
| 输入 | `sources/sympy/sympy/**/tests/test_*.py` |
| 解析 | AST 找 `assert` / `assertEqual` 中含 `limit`/`integrate`/`diff`/`solve` 的调用 |
| 输出 | `fixtures/sympy/{limits,integrals,...}.draft.json` |
| `upstream_input` | SymPy 源码片段 |
| `line` | 手工或规则转 giac；转不了则 `verify: parse_only` |
| CLI | `--module series.tests.test_limits` 等 |

---

## 7. `extract_rubi_problems.py`（规格摘要，待实现）

| 项 | 约定 |
|----|------|
| 输入 | Rubi 下载的纯文本 / Maxima 语法包 |
| 格式 | 每行或每条：`(integrand, var)` |
| 输出 | `fixtures/rubi/samples/<category>.draft.json` |
| `verify` | 固定 `integrate_derivative` |
| 抽样 | `--max N` · `--category 4.3.11` · `--manifest rubi/manifests/catalog.yaml` |

---

## 8. CI 集成

### 8.1 PR 门禁（默认）

```bash
cd giac-rs
cargo test --workspace
cargo test -p giac-conformance --test external_maxima enabled
# 仅跑 enabled: true 的 maxima fixture
```

### 8.2 可选 / nightly

```bash
cargo test -p giac-conformance --test external_rubi -- --ignored
cargo test -p giac-conformance --test external_sympy
```

### 8.3 与 giac-calculus 内嵌测试

`giac-calculus::limit::tests::maxima_rtest`（14 条）为 **crate 内联** 快测；与 `fixtures/maxima/limit.json` **应对齐**，但不必重复 id。长期：fixture 为唯一源，crate 测试 `include_str!` 或 build 脚本生成。

---

## 9. 许可与合规

| 上游 | 许可 | giac-rs 使用 |
|------|------|--------------|
| Maxima rtest | GPL-2.0+ | 抽取表达式进 JSON；注明 `source` |
| SymPy tests | BSD-3-Clause | 同上 |
| Rubi 题干 | MIT | 同上 |
| giac bin/check | GPL-3.0 | T5 主规格 |

**不复制：** Mathematica/Maple 专有测试框架；仅可复用 **数学表达式**（facts）。

---

## 10. 维护 checklist

1. 上游版本 bump → 重跑 extract → diff `*.draft.json` → 审核
2. 新 `enabled: true` → 必须绑定 `giac_issue` 或 phase 里程碑
3. 关闭 issue → 对应 fixture 行改 `enabled: true` 并 CI 绿
4. 失败时优先查 T1（SymPy），再 `assert_equiv`，最后查 [known-divergences.md](known-divergences.md)

---

## 参考

- [conformance-testing.md §3](conformance-testing.md) — `assert_equiv`
- [phase4-issues.md §2.6](phase4-issues.md) — SymPy 判定模式
- [GIAC-limit-maxima-upstream-alignment.md](issues/GIAC-limit-maxima-upstream-alignment.md) — 14 条 limit 子集
- [12000.org CAS integration tests](https://www.12000.org/my_notes/CAS_integration_tests/) — Rubi 横向基准
