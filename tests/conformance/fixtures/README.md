# Conformance fixtures

外部与内部 curated 测试用 JSON。规格见 [`.doc/external-test-resources.md`](../../../../.doc/external-test-resources.md)。

| 目录 | 来源 | 说明 |
|------|------|------|
| `_schema/` | — | `fixture-v1.schema.json` |
| `internal/` | giac bin/check | 索引文档，不复制 golden |
| `maxima/` | Maxima rtest | `extract_maxima_rtest.py` 输出 |
| `sympy/` | SymPy tests | `extract_sympy_tests.py`（待实现） |
| `rubi/` | Rubi 积分集 | 抽样子集，默认不进 PR 门禁 |

**工作流：** `*.draft.json`（全 disabled）→ 人工审核 → `*.json` → 部分 `enabled: true`。
