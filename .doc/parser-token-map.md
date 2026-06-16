# 解析器 Token / 优先级 / 多模式 / Context 映射

对照 giac-1.5.0：`input_lexer.ll` → `input_lexer.cc`（flex 2.5）、`input_parser.yy` → `input_parser.cc`（bison）。

Rust 目标 crate：`giac-parse`；MVP 仅实现 **Xcas 脚本子集**（`xcas_mode(0)`，无 TI/Maple/RPN 方言）。

---

## 1. 构建链（flex / bison）

| 步骤 | 命令 / 文件 | 产物 |
|------|-------------|------|
| 词法 | `flex input_lexer.ll` | `input_lexer.cc`, `lexer.h` |
| 语法 | `bison -p giac_yy -y -d input_parser.yy` | `input_parser.cc`, `input_parser.h` |
| 编译 | 与 `gen.cc` 等链入 `giac` 可执行文件 | — |

**CMake 现状：** 仓库内 **已提交** 生成后的 `input_lexer.cc` / `input_parser.cc`，日常构建 **不** 调用 flex/bison。修改 `.ll`/`.yy` 需本地安装 flex 2.5 + bison 后手动再生成。

**Rust 侧：** 不使用 flex/bison；用 `logos` + `lalrpop` 手写等价 grammar，token 名以下表 `Rust Token` 列为准。

---

## 2. 完整 Token 清单

### 2.1 字面量与标识符

| C++ Token | 典型词法 | 语义 | Rust Token |
|-----------|----------|------|------------|
| `T_NUMBER` | `123`, `3.14` | 整数/浮点字面量 | `Number` |
| `T_SYMBOL` | `x`, `foo` | 标识符 / 函数名 | `Ident` |
| `T_LITERAL` | `pi`, `∞`, `undef` | 预定义常量 | `Literal` |
| `T_DIGITS` | `0x1A` 等 | 进制数字 | `Digits` |
| `T_STRING` | `"hello"` | 字符串 | `Str` |
| `T_EXPRESSION` | 特殊表达式 token | 预解析表达式 | — |
| `T_END_INPUT` | EOF | 输入结束 | `Eof` |

### 2.2 运算符与标点

| C++ Token | 字符/形式 | 构造 AST | Rust Token |
|-----------|-----------|----------|------------|
| `T_PLUS` | `+` | `at_plus` | `Plus` |
| `T_MOINS` | `-`（二元） | `symb_plus(a, -b)` | `Minus` |
| `T_MOINS38` | `-`（TI-38 二元） | 同上 | `Minus38` |
| `T_FOIS` | `*` | `at_prod` | `Star` |
| `T_IMPMULT` | （无字符，语法优先级） | 隐式乘 `2x` | （parser 规则） |
| `T_DIV` | `/` | `at_div` | `Slash` |
| `T_MOD` | `mod` / `%` | `normalmod` / `at_mod` | `Mod` |
| `T_POW` | `^` | `at_pow` | `Caret` |
| `T_SQ` | `**` | `symb_pow` | `StarStar` |
| `T_PRIME` | `'` | `at_derive` | `Prime` |
| `T_FACTORIAL` | `!` | `at_factorial` | `Bang` |
| `T_NOT` | `not` | 一元 not | `Not` |
| `T_TEST_EQUAL` | `==`, `<=`, `>=`, `<`, `>` | 关系运算符 | `RelOp` |
| `T_EQUAL` | `=` | 方程 `at_equal` | `Eq` |
| `T_AND_OP` | `and` / `&&` | 逻辑与 | `And` |
| `T_VIRGULE` | `,` | 序列/参数分隔 | `Comma` |
| `T_SEMI` | `;` | 语句结束 | `Semi` |
| `T_AFFECT` | `:=` | `symb_sto` | `Assign` |
| `T_MAPSTO` | `->` | 程序映射 | `Mapsto` |
| `T_DEUXPOINTS` | `:` | 键值/比例 | `Colon` |
| `T_DOUBLE_DEUX_POINTS` | `::` | 特殊分隔 | `ColonColon` |
| `T_QUOTE` | `'` / `"` 引号控制 | `symb_quote` | `Quote` |
| `T_COMPOSE` | `@` / `o` | 复合 | `Compose` |
| `T_DOLLAR` | `$` | 电子表格 | `Dollar` |
| `T_DOLLAR_MAPLE` | `$`（Maple） | `symb_dollar` | —（不做） |
| `T_PIPE` | `\|` | 条件概率等 | `Pipe` |
| `T_INTERROGATION` | `?` | 帮助/索引 | `Question` |

### 2.3 容器与下标

| C++ Token | 形式 | subtype / 说明 | Rust |
|-----------|------|----------------|------|
| `T_BEGIN_PAR` / `T_END_PAR` | `( )` | 分组 / 函数参数 | `LParen` / `RParen` |
| `T_VECT_BEGIN` | `[` | 向量开始 | `LBracket` |
| `T_VECT_END` | `]` | 向量结束 | `RBracket` |
| `T_VECT_DISPATCH` | `[` 后关键字 | 见 §2.5 | `ContainerKind` |
| `T_SET_BEGIN` / `T_SET_END` | `{ }` | 集合 | `LBrace` / `RBrace` |
| `T_INDEX_BEGIN` | `[` 后缀下标 | `symb_at` | `LBracket` |
| `T_MATRICE_BEGIN` / `T_MATRICE_END` | `[[` `]]` | 矩阵 | `LMat` / `RMat` |
| `T_ROOTOF_BEGIN` / `T_ROOTOF_END` | `rootof[...]` | 代数扩域 | `RootOf` |
| `T_SPOLY1_*` / `T_POLY1_*` | 稀疏/稠密多项式 | 内部多项式 | 后置 |
| `T_ASSUME_BEGIN` / `T_ASSUME_END` | 假设块 | `_ASSUME__VECT` | `Assume` |

### 2.4 集合运算

| Token | 运算 |
|-------|------|
| `T_UNION` | 并集 |
| `T_INTERSECT` | 交集 |
| `T_MINUS` | 集合差 |
| `T_INTERVAL` | 区间 `a..b` |

### 2.5 `T_VECT_DISPATCH` 子类型（词法设置 subtype）

| subtype 常量 | 语法示例 | Rust `Expr` 变体 |
|--------------|----------|------------------|
| `_SEQ__VECT` | `[a,b,c]` 默认 | `Seq` |
| `_LIST__VECT` | `list[a,b]` / `{a,b}` 模式 | `List` |
| `_SET__VECT` | `{a,b}` | `Set`（后置） |
| `_MATRIX__VECT` | `[[1,2]]` | `Matrix` |
| `_POLY1__VECT` | 多项式系数 | 内部 `Poly1` |
| `_ASSUME__VECT` | 假设 | `Context.assumptions` |
| `_PNT__VECT`, `_CURVE__VECT`, … | 几何对象 | 后置 / geo |

### 2.6 控制流与程序（MVP 子集加粗）

| Token | 关键字 | MVP |
|-------|--------|-----|
| `T_FOR`, `T_FROM`, `T_TO`, `T_DO`, `T_BY`, `T_IN` | `for` | **✓** `test_prog` |
| `T_WHILE`, `T_MUPMAP_WHILE` | `while` | **✓** |
| `T_IF`, `T_THEN`, `T_ELSE`, `T_ELIF`, `T_IFTE` | 条件 | △ flanex |
| `T_REPEAT`, `T_UNTIL` | `repeat` | ✗ |
| `T_RETURN`, `T_BREAK`, `T_CONTINUE` | 控制 | △ |
| `T_TRY`, `T_CATCH`, `T_TRY_CATCH`, `T_IFERR` | 异常 | ✗ |
| `T_LOCAL`, `T_PROC`, `T_PROGRAM`, `T_BLOC` | 程序块 | **✓** 部分 |
| `T_SWITCH`, `T_CASE`, `T_DEFAULT` | switch | ✗ |

### 2.7 TI / Maple / RPN 专用（MVP 不做）

| Token 组 | 模式 | 说明 |
|----------|------|------|
| `TI_STO`, `TI_SEMI`, `TI_DEUXPOINTS`, `TI_FOR`, `TI_WHILE`, … | `xcas_mode==3` / TI | 存储 `→`、分号语义不同 |
| `T_RPN_*`, `T_STACK`, `T_ACCENTGRAVE` | `rpn_mode` | 逆波兰 |
| `T_MAPLELIB`, `T_DOLLAR_MAPLE` | Maple | `$` 求和等 |
| `T_CASE38`, `T_MOINS38`, `T_NEG38`, `T_UNARY_OP_38` | `abs_calc_mode==38` | NumWorks/TI-38 |
| `T_LOGO` | Logo 方言 |  turtle |

### 2.8 一元函数 token

| Token | 说明 |
|-------|------|
| `T_UNARY_OP` | 词法表查到的函数名 → `gen(at_xxx, arity)` |
| `T_UNARY_OP_38` | TI-38 一元 |
| `T_TYPE_ID` | 类型名 `matrix`, `list` 等 |
| `T_UNIT` | 单位 `5_m` |

---

## 3. 运算符优先级（bison，低 → 高）

行号对应 `input_parser.yy` 声明顺序：**越靠后优先级越高**。

| 层级 | Token / 组 | 结合性 | 示例 |
|------|------------|--------|------|
| 1 | `T_AFFECT` (`:=`) | 右 | `x:=1` |
| 2 | `TI_SEMI` | 左 | TI 语句 |
| 3 | `T_VIRGULE` | 左 | `a,b` |
| 4 | `T_AND_OP` | 左 | `a and b` |
| 5 | `T_EQUAL` | 右 | 方程 |
| 6 | `T_TEST_EQUAL` | 左 | `a==b` |
| 7 | `T_UNION` | 左 | 集合并 |
| 8 | `T_INTERSECT` | 左 | 交集 |
| 9 | `T_INTERVAL` | 左 | `1..10` |
| 10 | `T_PLUS`, `T_MOINS`, `T_MOINS38` | 左 | 加减 |
| 11 | `T_FOIS`, **`T_IMPMULT`** | 左 | `*`, 隐式乘 |
| 12 | `T_DIV` | 左 | `/` |
| 13 | `T_MOD` | 非结合 | `mod` |
| 14 | `T_POW` | 右 | `^` |
| 15 | `T_FACTORIAL` | 非结合 | `!` |
| 16 | `T_SQ` | 左 | `**` |
| 17 | `T_UNARY_OP`, `T_NEG38`, `T_NOT` | 非结合 | 一元 |
| 18 | `T_COMPOSE` | 左 | 复合 |

**隐式乘法规则（`T_IMPMULT`）：**

```text
T_NUMBER symbol_or_literal           →  2*x  或  2*(x+1) 经语法树
T_NUMBER symbol_or_literal T_POW ... →  2*x^3
T_NUMBER T_UNARY_OP ( exp )          →  2*sin(x)
```

测试触发：`check/testcas`、`bin/test_cas_basic` 中系数紧邻变量。

---

## 4. 多计算器模式

模式标志存于 `context` / `globalptr`（`global.h`），词法通过 `yyextra` 读取。

| 标志 |  accessor | 值 | 词法/语法影响 | MVP |
|------|-----------|-----|---------------|-----|
| **xcas_mode** | `xcas_mode(ctx)` | 0=Xcas, 1=Maple, 2=Mupad, 3=TI | `;`/`:`  token 类型；`i` 是否 sqrt(-1) | **0** |
| **calc_mode** | `calc_mode(ctx)` | 0/1/38… | 区间、矩阵 1-based 等 | **0** |
| **abs_calc_mode** | `abs_calc_mode(ctx)` | 38=NumWorks | 专用 `-38` token | 不做 |
| **rpn_mode** | `rpn_mode(ctx)` | 0/1 | RPN 栈语法 | 不做 |
| **python_compat** | `python_compat(ctx)` | bool | `()` → python_list | 不做 |
| **i_sqrt_minus1** | `i_sqrt_minus1(ctx)` | 0/1 | 解析后二次扫描 | **跟随测试** |

**测试脚本默认初始化（check 套件首行）：**

```giac
cas_setup(0,0,0,1,0,[1e-10,1e-17],12,[1,50,0,25],0,0,0),xcas_mode(0);
```

含义见 §5。Rust conformance harness 应在跑 golden 前注入等价 `Context` 初始化。

---

## 5. `cas_setup` → `Context` 映射

实现：`prog.cc::cas_setup(const vecteur & v, GIAC_CONTEXT)`（约 6654 行）。

参数 **至少 7 项**；测试与 segfault 用例使用 **12 项**。

| 索引 | giac 参数 | 类型 | Context 字段（Rust） | 测试典型值 |
|------|-----------|------|----------------------|------------|
| v[0] | approx_mode | bool/int | `approx_mode: bool` | `0` 精确 |
| v[1] | complex_variables | bool | `complex_variables: bool` | `0` |
| v[2] | complex_mode | bool | `complex_mode: bool` | `0` |
| v[3] | angle_mode | 0=rad, 1=deg, 2=grad | `angle_mode: AngleUnit` | `1` → **degrees** |
| v[4] | format | int | `scientific_format`, `integer_format` | `0` |
| v[5] | epsilon | float 或 `[eps, proba_eps]` | `epsilon: f64`, `proba_epsilon: f64` | `[1e-10, 1e-17]` |
| v[6] | decimal_digits | int | `float_digits: u32` | `12` |
| v[7] | `[threads, MAX_REC, debug, eval_level]` | vecteur | `threads`, `max_recursion`, `debug`, `eval_level` | `[1,50,0,25]` |
| v[8] | increasing_power | bool | `increasing_power: bool` | `0` |
| v[9] | withsqrt | bool | `withsqrt: bool` | `0` |
| v[10] | all_trig_sol | bool | `all_trig_sol: bool` | `0` |
| v[11] | integer_mode | bool | `integer_mode: bool` | `0` |

**Rust 初始化示例：**

```rust
Context {
    approx_mode: false,
    complex_variables: false,
    complex_mode: false,
    angle_mode: AngleUnit::Degrees,
    epsilon: 1e-10,
    proba_epsilon: 1e-17,
    float_digits: 12,
    eval_level: 25,
    max_recursion: 50,
    increasing_power: false,
    withsqrt: false,
    all_trig_sol: false,
    integer_mode: false,
    xcas_mode: 0,
    ..Default::default()
}
```

---

## 6. `assume` / `purge` → `Context`

### 6.1 `assume`（`usual.cc::giac_assume`）

| 形式 | 行为 | Rust |
|------|------|------|
| `assume(x, real)` | `sto(_ASSUME__VECT, x)` 类型实数 | `assumptions.push(Real(x))` |
| `assume(x, complex)` | 复数假设 | `Complex(x)` |
| `assume(x > 0)` | 符号假设 `assumesymbolic` | `Relation(x, GT, 0)` |
| `assume([...])` | 合取 | 多条 assumption |
| `additionally` | 追加假设 | `add_assumption` |

存储：变量绑定为 `_ASSUME__VECT` 子类型向量，化简/积分时读取（`sym2poly.cc::check_assume`）。

### 6.2 `purge`（`rpn.cc::_purge`）

| 形式 | 行为 | Rust |
|------|------|------|
| `purge(x)` | 删除变量及其假设 | `vars.remove(x)` + 清 assumption |
| `purge()` | 清全部用户变量 | `vars.clear()` |

测试来源：`check/flanex`, `check/testgeo`。

---

## 7. Rust `giac-parse` MVP 对照摘要

| 能力 | 必须 | 可推迟 |
|------|------|--------|
| Token 子集 §2.1–2.3 | ✓ | TI/RPN/Maple §2.7 |
| 优先级 §3 | ✓ | Logo/spreadsheet |
| `cas_setup` 12 参数 | ✓ harness 注入 | GUI 动态修改 |
| `assume`/`purge` | ✓ | `additionally` 全分支 |
| 隐式乘 | ✓ | — |
| `for`/`while` | ✓ | `switch`/`try` |
