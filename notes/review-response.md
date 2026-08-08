# looplook — Bioconductor 审稿回复（round 2, Bisaloo, 2026-08-05）

## 对应审稿人原话

> Is there no existing function/package to read `.bed` / `.bedpe` files already in Bioconductor?

## 可直接粘贴的回复（英文）

Thanks for the question about existing `.bed`/`.bedpe` readers. We evaluated
`rtracklayer::import()`, which does support both formats (BED → `GRanges`; BEDPE →
`Pairs`). Your question also prompted us to check that our own consolidated output
is standard-compliant, and we fixed an inconsistency it exposed.

**Standard-compliant BEDPE output.** `consolidate_chromatin_loops()` now writes
standard 10-column BEDPE (columns 9–10 are `strand1`/`strand2`, set to `"."` as the
data are not stranded), with the documented consensus metadata `n_members`/`n_reps`
appended as extension columns 11–12. This makes the exported files directly
readable by `rtracklayer::import(..., format = "BEDPE")` (verified: returns a
`Pairs` object with the score and the extension columns preserved) as well as by
standard BEDPE tooling.

**Why a thin custom reader is still retained.**

- **Score semantics.** `bedpe_to_gi()` auto-detects the score column (columns 7–8,
  by numeric ratio) and warns when the detected column resembles p-values.
  rtracklayer's BEDPE import always treats column 8 as a fixed score, with no such
  detection — non-standard user layouts would be silently mis-read.
- **Extension columns.** `bedpe_to_gi()` recovers the documented extension
  columns 11–12 as `n_members`/`n_reps` (integer-gated, so arbitrary 12-column
  files are not mislabelled), making the consolidated-output round-trip lossless.
  rtracklayer preserves those columns only as unnamed metadata (`NA.`/`NA..1`).
- **Heterogeneous ecosystem input.** `rtracklayer::import(..., format = "BEDPE")`
  is scoped to strict, headerless 10-column BEDPE. Real published loop files are
  not always in that layout. For example, the K562 ChIA-PET loop calls shipped
  with the `GenomicInteractions` package carry a header row and tool-specific
  column names; rtracklayer rejects the file (verified error: `scan() expected an
  'integer', got 'start_left'`), whereas `bedpe_to_gi()` skips the header by
  explicit header tokens and reads the file (verified: 64,565 loops). We do not
  claim universal coverage: files whose coordinate columns are not in positions
  1–6 (e.g. the Seitan et al. 2013 interaction file, which carries an ID prefix)
  remain outside the current reader's scope and would require column remapping;
  they are documented as such. To be clear, this is not a criticism of rtracklayer
  — those files fall outside its documented BEDPE scope. It simply reflects that a
  lenient reader (skip headers, auto-detect the score column, tolerate variable
  column counts) has practical value for the loop-file ecosystem, in addition to
  the standard files that rtracklayer already handles.
- **Performance.** The parser is `data.table::fread`-based; on a representative
  genome-wide BEDPE (~500k rows, ~30 MB) file reading is roughly an order of
  magnitude faster than rtracklayer's object-building import, which matters for
  million-row whole-genome HiChIP inputs.

The parsing layer is deliberately thin, is covered by the `data_processing` test
suite (45 test blocks / 93 assertions), and its public boundary is the standard
Bioconductor `GInteractions` class — the file-I/O step is internal, not a
re-implementation of ecosystem object models.

For plain 3-column BED (blacklist, ROI, peaks), `rtracklayer::import("BED")`
produces field-identical output to our reader. Our reader additionally tolerates
header rows and emits narrow/wide-interval diagnostics; switching the BED path to
rtracklayer would drop those two behaviors and require updating the corresponding
tests. We are happy to make that switch if you prefer it, at the cost of those
behaviors — please let us know.

## 论证的实测证据（回复中引用的每条都有据可查）

### BEDPE round-trip（rtracklayer 1.70.1 实测）

```r
# 当前 looplook 导出（标准 10 列前缀 + 扩展列 11/12）
chr1	0	100	chr1	200	300	.	42	.	.	2	3
chr1	500	600	chr1	800	900	.	7	.	.	1	2

rtracklayer::import(file, format = "BEDPE")
# OK: Pairs 对象，score 正确，扩展列保留为 mcols "NA."/"NA..1"
```

改动前（10 列、9/10=n_members/n_reps）rtracklayer 报错
`strand values must be in '+' '-' '*'`——审稿人质疑"为何非标准"成立，已修复。

修复：`.consolidate_export()`（`R/data_processing.R:1248`）现导出 12 列：
chr1, start1, end1, chr2, start2, end2, name, score, **strand1=".", strand2="."**,
n_members, n_reps。数值逐列核对与 `GInteractions` 对象一致（start 0-based、
score、n_members col11、n_reps col12）。

`n_members`/`n_reps` 为文档化元数据：`@return`（`R/data_processing.R:664-669`）
与 vignette 已同步为"标准 10 列 + 扩展列 11/12"。

### 真实发布数据的异质性（rtracklayer 1.70.1 / R 4.5.1 实测）

rtracklayer 的 BEDPE 解析器面向**严格、无表头的 10 列**布局（9/10 列必须是
strand）。但真实发布的 Hi-C/ChIA-PET 循环文件经常不遵循该布局。用
`GenomicInteractions` 包（已发布数据的参考实现）自带的两份真实文件实测：

| 文件 | 来源 | 布局 | rtracklayer `import("BEDPE")` | looplook `bedpe_to_gi` |
|---|---|---|---|---|
| `k562.rep1.cluster.pet3+.txt` | K562 ChIA-PET 循环（Seitan 等） | 9 列 + 表头，坐标在 1-6 | **报错** `got 'start_left'` | **OK，64,565 loops**（显式 token 跳表头） |
| `Seitan2013_WT_100kb_interactions.txt` | Seitan 2013 发表的 100kb 互作 | 20 列 + 表头，坐标在 3-6 | **报错** `got 'PeakID.1.'` | **不支持**（坐标不在 1-6，需列重排） |

**严谨的表述边界**：
- K562 属"表头 + 工具特有列名、坐标在 1-6"——looplook 能读、rtracklayer 不能，
  这是**真实的、可复核的功能差异**。
- Seitan 属"坐标带 ID 前缀、需列重排"——超出当前 reader 范围，**如实写明"按需
  可扩展、暂不支持"**，不假装能读。
- 两份文件均非标准 10 列 BEDPE，落在 rtracklayer 文档化 BEDPE 范围之外——
  rtracklayer 拒绝它们是其**符合文档化范围**的行为，不是 rtracklayer 的缺陷。
  只论证"宽容解析器（跳表头 + 分数自动检测 + 灵活列数）有实用价值"，
  不做"rtracklayer 不好用"的过度推论。

配套实现（`R/data_processing.R` `bedpe_to_gi`）：
- 表头检测仅凭**显式 token**（`chrom_left/start_left/...`），不因数值失败跳过
  （保留 "malformed first row must error" 语义，与 `read_simple_bed` 一致）。
- 扩展列 11/12 恢复按**整数化门控**（`.numeric_ratio ≥ 0.5` 且 `== floor`），
  避免任意 12 列文件被误标。

### 分数自动检测 / p-value 告警（代码 + 测试）

- `bedpe_to_gi()`：`.numeric_ratio(df[[col]]) >= 0.5` 自动判定第 7/8 列；`.is_pvalue_like()` 触发告警
  （`R/data_processing.R`）。
- 测试：`tests/testthat/test-data_processing.R:475` "warns when auto-detected score looks like p-values"。
- rtracklayer BEDPE 固定把第 8 列当 score，无检测。

### 性能（实测，500k 行 / 30.5 MB）

| 解析方式 | 耗时 |
|---|---|
| `data.table::fread` | 0.05 s |
| `rtracklayer::import(format="BEDPE")` | 0.80 s（约 16×） |

### BED 路径对比（实测，`example_peaks.bed` 300 区间）

| 场景 | looplook `read_simple_bed` | rtracklayer `import("BED")` |
|---|---|---|
| 纯 3 列 BED | 300 ranges | 300 ranges，seqnames/starts/ends **逐一 identical** |
| 带表头行（start/end） | 检测并跳过 → 正常读 | **报错** `scan() expected 'an integer', got 'start'` |
| 窄区间 <10bp | 读出 + 宽度告警 | 读出，无告警 |
| start > end | 报错 "start >= end" | 报错（消息不同） |
| 4 列（name） | 丢弃 name | 保留为 mcols |
| 空/纯空白文件 | 报错 | 返回 0 ranges |

含义：BED 让步**不影响标准文件结果**，但会丢失表头容忍与宽度诊断两个行为，并需改写
`tests/testthat/test-coverage-boost.R:628,639,647` 等断言具体消息的测试。

### 测试基数

`tests/testthat/test-data_processing.R`：45 个 `test_that` 块 / 93 条 `expect_*` 断言
（实测分布：expect_equal 37, expect_false 10, expect_error 9, expect_s4_class 8,
expect_true 8, expect_null 3, expect_warning 3, expect_setequal 2,
expect_no_message 1, expect_no_warning 1；无注释/字符串假阳性）。

## 相对初稿的措辞改动

1. **论证方向反转**：审稿人质疑"为何非标准"成立 → 我们把导出改为标准 10 列
   BEDPE + 扩展列 11/12。回复从"守住非标准格式"（自曝弱点）改为"我们修正了
   格式，rtracklayer 现在能直接读"（进攻性陈述）。round-trip 从"报错"变"成功"。
2. 自定义 reader 的理由从 3 条扩为 4 条：分数自动检测/告警 + 扩展列恢复
   （rtracklayer 只保留为无名 `NA.`/`NA..1`）+ **生态异质性** + 性能。
3. 性能用 "roughly an order of magnitude"，不写死秒数。
4. 测试数用实测的 "45 test blocks / 93 assertions"。
5. 补 "not a re-implementation of ecosystem object models"，正面回应复用纪律。
6. **两条声称都做成真的**（初稿曾两头空，已被内部审查发现）：
   - "扩展列恢复"：初稿声称 `bedpe_to_gi` 恢复 n_members/n_reps，实际没实现——
     **已实现**（整数化门控，往返无损，`test-data_processing` +2 测试锁定）。
   - "生态异质性"：初稿举例的两份真实文件 looplook 自己读不了——**K562 已实现
     表头容忍**（显式 token，与 `read_simple_bed` 同模式），现在 looplook 能读、
     rtracklayer 不能，成为真实可复核的差异；**Seitan 如实写明"坐标在 3-6、需列
     重排、暂不支持"**，不假装能读。审稿人若喂文件给 `bedpe_to_gi`，K562 能过、
     Seitan 得到明确说明，不会再被"两头空"反将。
7. BED 从"无条件让"改为"field-identical 但会丢表头容忍 + 宽度诊断 + 需改测试"，
   姿态配合但量化代价，交由审稿人决定。

## 战略定位

- **导出格式：修正为标准 BEDPE**（10 列前缀 + 扩展列 11/12）——把审稿人可能
  质疑的"非标准格式"弱点直接修掉，rtracklayer 现可读回。
- **BEDPE 解析器：守住**（分数自动检测/告警 + 扩展列无损往返 + K562 式表头
  容忍 + 性能），且**声称与实际能力一致**——每一条都有实测和测试兜底。
- **BED 解析器：可让，但有明确代价**（表头容忍、宽度诊断、测试重写）。
- **回复姿态**：不否认现成函数，不宣称"更好"，只陈述"评估过 + 实测证据 +
  具体代价"，且主动把发现的不一致修掉——这是审稿人最想看到的回应。

---

# Part 2: `"Set 2"` / RColorBrewer 依赖（审稿人 A 表 #2）

## 对应审稿人原话

> Note that `"Set 2"` is provided by `grDevices::palette()` so the RColorBrewer
> dependency might not be necessary.
> （引文指向 `R/analysis.R` 的 Set2 硬编码处）

## 可直接粘贴的回复（英文）

Thanks — you are right that the qualitative ColorBrewer family (including
`"Set 2"`, `"Dark 2"`, `"Paired"`) is available from `grDevices::palette.colors()`
since R ≥ 4.0. We still keep `RColorBrewer` in `Imports`, for three reasons:

1. **Diverging/sequential palettes are not in grDevices.**
   `grDevices::palette.pals()` covers only the 8 qualitative palettes; it has no
   diverging family such as `"PuOr"`. The package uses `"PuOr"` for the
   log2-fold-change colour scale in the GO cnet plot and the PPI network
   (`ggplot2::scale_color_distiller(palette = "PuOr")` in `run_go_enrichment()`
   and `run_ppi_analysis()`), which grDevices cannot provide.

2. **The dependency is present transitively regardless.**
   `scale_color_distiller()` is rendered through the `scales` package, whose
   `Imports` include `RColorBrewer`. Removing our direct use would not remove
   `RColorBrewer` from the dependency graph — it would only turn an explicit
   requirement into an implicit one.

3. **The public `color_palette` API documents RColorBrewer names.**
   `get_colors()` resolves arbitrary documented palette names (e.g. `"Set2"`,
   `"Dark2"`, `"Paired"`) via `RColorBrewer::brewer.pal.info`, and the
   hard-coded qualitative calls use the same ColorBrewer family. Routing some
   calls through `grDevices` and others through RColorBrewer would split the
   single source of truth without saving a dependency.

Declaring `RColorBrewer` explicitly is therefore clearer than relying on
`scales`' internals.

## 论证的实测证据

### grDevices 自带范围（R 4.4.1 实测）

`grDevices::palette.pals()` 返回 16 个定性调色板：R3, R4, ggplot2, Okabe-Ito,
Accent, **Dark 2, Paired, Pastel 1, Pastel 2, Set 1, Set 2, Set 3**, Tableau 10,
Classic Tableau, Polychrome 36, Alphabet。

- `"Set 2"`/`"Dark 2"`/`"Paired"`：grDevices 有 ✅
- **`"PuOr"`：grDevices 无 ❌**（发散型不在内置定性列表中）

### PuOr 实际用法（`R/analysis.R`）

- `analysis.R:1418`：`run_go_enrichment()` 的 GO cnet 图
  `scale_color_distiller(palette = "PuOr", name = "Log2FC")`
- `analysis.R:1593`：`run_ppi_analysis()` 的 PPI 网络
  `scale_color_distiller(palette = "PuOr", name = "LFC")`

### 依赖链（R 4.4.1 实测）

| 包 | Imports 是否含 RColorBrewer |
|---|---|
| `ggplot2` | 否 |
| `scales`（ggplot2 依赖） | **是** |

`scale_color_distiller()` → `scales::brewer_pal()` → `RColorBrewer::brewer.pal()`。
即使用户不用定性色板，`"PuOr"` 这一处就保证 RColorBrewer 必然进入依赖图。

### 包内 RColorBrewer 用法汇总

- 定性硬编码：`analysis.R:1639`（Dark2）、`analysis.R:1728,1730`（Set2）
- 通用解析：`utils.R:1926-1929` `get_colors()` 经 `brewer.pal.info` 接受任意
  RColorBrewer 调色板名；`utils.R:1917-1924` 注释已写明保留理由（a/b/c 三条）
- 公开 API：`color_palette` 默认 `"Set2"`（`analysis.R:2912`、`annotation.R:928`）

## 相对 Part 1 的措辞注意点

1. 首句**先承认审稿人正确**（Set 2 确实在 grDevices），不反驳事实，只补全
   "但依赖去不掉"的完整理由。
2. 理由 1（PuOr）是**决定性**的：grDevices 无发散型色板，而包内两处用到。
3. 理由 2（scales 传递依赖）封住"那你干嘛显式声明"的追问：显式声明比隐式
   依赖更清晰，且删了也省不掉。
4. 理由 3（公开 API 以 RColorBrewer 命名为准）说明不切换的原因——切了会
   分裂唯一事实源，且不省任何依赖。

## 战略定位

- **RColorBrewer：守住，理由充分**（PuOr 发散型 grDevices 无 + scales 传递依赖
  删不掉 + 公开 API 命名锚定）。
- **姿态**：不写 "you're wrong"，写 "you're right, and here is why the dependency
  stays anyway"——先给审稿人认可，再补全上下文。

---

# Part 3: 总建议（可读性 / 可维护性 / 抽 utils / 去重复）

## 对应审稿人原话

> Beyond the specific comments in the review below, I would strongly encourage
> you to make any edits possible to improve the package readability and
> maintainability: define accessory utils functions for commonly repeated tasks
> and use these functions wherever possible, refactor duplicated code, etc.

## 可直接粘贴的回复（英文）

Thanks for this overall suggestion. In this round we made a targeted pass over
the largest duplication in the codebase and further consolidated shared helpers:

- **Workbook-export boilerplate.** The three Excel exporters
  (`.export_annotation_excel()`, `.refine_export_workbook()`,
  `.chromatin_export_excel()`) previously repeated the same
  `createWorkbook → addWorksheet/writeData → saveWorkbook+tryCatch` skeleton.
  We extracted two shared utilities, `.add_sheet()` (null/empty handling,
  optional column stripping, worksheet creation) and `.save_workbook()`
  (path construction + failure fallback), removing ~80 lines of near-identical
  code. Sheet names, sheet order, column stripping, and warning messages are
  byte-identical to the previous output, and the workbook-structure tests pass
  unchanged.
- **Output-directory handling.** The `write_output = TRUE` → validate/create
  `out_dir` contract was duplicated across the pipeline functions; it is now a
  single `.ensure_out_dir()`, which also fixed an inconsistency where one
  refinement path validated but did not create the directory.
- **Scalar argument validation.** Manual inline guards for `seed`, `score_col`
  and `min_score` were replaced by the existing `.assert_scalar_count()` /
  `.assert_scalar_number()` helpers (this also removed a case where `min_score`
  was validated by hand in one function and by the helper in another). The
  remaining manual checks were intentionally left as-is where a helper cannot
  express their semantics exactly (e.g. a strictly-positive bound, or a
  quantile-specific message) — see below.

We also rely on existing helper discipline throughout the package: shared
splitting/parsing and annotation helpers are reused broadly (e.g. `extract_genes`
in 41 call sites, `.with_known_upstream_noise_suppressed` in 45, a common
`looplook_metadata` builder, and a consistent `.assert_*` validation family),
with the package kept as small, single-purpose functions (210 across the R/ tree).
All of these changes are covered by the existing test suite (773 assertions across
the files exercised in this round, 0 failures), so the refactoring is
behaviour-preserving.

## 实测依据

### Excel 样板重构的量化

| 函数 | 重构前 | 重构后 |
|---|---|---|
| `.export_annotation_excel` | 35 行 | 12 行 |
| `.refine_export_workbook` | 45 行 | 15 行 |
| `.chromatin_export_excel` | 40 行 | 12 行 |

净减约 80 行；三个函数共享 `.add_sheet()` / `.save_workbook()`（`R/utils.R`）。

### 行为保持的证据

- sheet 名、顺序、列清洗（a1_id/a2_id、loop_genes、SANKEY_RAW_GENES 等）、
  文件名后缀（`_Basic_Results.xlsx` / `_Refined_Results.xlsx` /
  `_Chromatin_Results.xlsx`）、warning 前缀全部逐字保留。
- `require_rows = FALSE` 对应原 `!is.null()`；`require_rows = TRUE` 对应原
  `!is.null() && nrow() > 0`——空值语义逐一映射。
- 工作簿结构测试（`test-refine.R`、`test-review-fixes.R`）通过不变。

### 既有复用纪律（写进回复的数字）

- `extract_genes`：41 处调用（`R/annotation.R` 等）
- `.with_known_upstream_noise_suppressed`：45 处调用（S4Vectors/ggplot2 已知
  上游噪音统一抑制）
- `.build_looplook_metadata`：6 处调用
- `.assert_*` 校验族（`.assert_scalar_count`、`.assert_nonempty_string`、
  `.assert_file_exists` 等）贯穿全部公共函数
- `.make_log_message(quiet)` 工厂：8 个 pipeline 函数共用
- 全包 210 个函数（analysis 31 / annotation 92 / utils 53 / visualization 14 /
  data_processing 20），小函数拆分良好

### 标量参数校验收敛到 `.assert_*`（实测依据）

迁移 3 处（干净 1:1，测试安全）：

| 参数 | 位置 | 迁移写法 | 备注 |
|---|---|---|---|
| `seed` | analysis.R:196 | `.assert_scalar_count(seed, "seed", min = 1)` | 消息仍含 "positive integer"，锁定测试（test-coverage-boost:1112）通过 |
| `score_col` | data_processing.R:175 | `.assert_scalar_count(score_col, "score_col", min = 1)` | 无测试锁定 |
| `min_score` | data_processing.R:364 | `.assert_scalar_number(min_score, "min_score")` | 顺带消掉"同参数一处手写、一处 .assert"的双处不一致 |

**刻意保留的 5 处**（helper 表达不了原语义，或消息被测试锁定，迁移会变消息/变行为）：

- `threshold`(quantile) annotation.R:3327：`[0,1]` 区间虽可用 `min=0,max=1`，但会丢
  "For quantile mode" 诊断上下文。
- `threshold`(绝对) annotation.R:3345：丢 "non-negative" 措辞。
- `ppi_species_id` analysis.R:1471：`.assert_scalar_count(min=1)` 更严——非整数当前
  静默 `as.integer()` 截断，迁移后直接报错（行为变更，虽更正确但非纯机械）。
- `bw_ratio_threshold` annotation.R:4616：要求严格 > 0，而 `.assert_scalar_number`
  的 `min` 是闭区间，`min=0` 会让 0 通过（行为变更）。

- `Expression threshold`（utils.R:1193 `.passes_expression_threshold`）：同族但消息
  "non-negative" 被 test-unit-utils.R:21 锁定——若迁到 `.assert_scalar_number`
  消息变 ">= 0" 即挂。这是"每处须逐个核对消息锁定"的实证。

### 测试基线

本轮涉及的全部相关测试文件：test-refine 30、test-review-fixes 248、
test-data_processing 93、test-regression 65、test-coverage-boost 295、
test-unit-utils 42，合计 **773 条断言，0 失败**。

## 相对初稿的措辞注意点

1. 不写"我们重写了整个包"——只陈述"对最大重复点做定向重构"，克制、可信。
2. "behaviour-preserving"要有硬证据：sheet 名/顺序/清洗/文件名逐字相同 +
   测试全绿，而不是空口承诺。
3. 复用纪律用具体数字（×41、×45、×6、210 函数），展示这不是口头承诺而是
   既有工程习惯。
4. 不承诺"未来继续重构"——避免给审稿人开空白支票；若被追问再说。

## 战略定位

- **已做**：三个定向重构——Excel 样板（`.add_sheet()`/`.save_workbook()`，
  ~80 行）、out_dir 校验（`.ensure_out_dir()`，顺带修一处不一致）、标量参数
  校验收敛（`.assert_*` 补齐 3 处）。全部行为逐字保持、测试锁定。
- **刻意不做**：低 ROI 项（空 GRanges/GInteractions 构造数处）；以及 helper 表达不了
  原语义的 5 处标量校验（threshold×2、ppi_species_id、bw_ratio_threshold、Expression
  threshold）——迁移会变消息/变行为。评审窗口期动低价值代码 = 回归风险换取几乎没有的收益。
  注：早期清单中的 `requireNamespace` 守卫与 `toupper` 重复调用，在本轮后续
  已分别收敛（见 Part 6 与前述 toupper 段落），不属"刻意不做"范围。
- **姿态**：展示"评估过 + 做了最高价值的 + 有测试兜底"，而不是把包重写一遍；
  刻意保留项在回复中主动说明理由（这就是审稿人想要的"评估过"证据）。

---

# Part 4: `trimws()` 标识符规范化（审稿人 A 表剩余风格项）

## 对应审稿人原话

> you keep using `trimws()` around column names (e.g. analysis.R L309). Do you
> have a valid reason to think there might be extra spaces here?

## 可直接粘贴的回复（英文）

On the `trimws()` calls around identifiers: yes, there is a concrete reason.
They are applied **at the input boundary** to user-supplied identifiers — gene
symbols, semicolon-split gene lists, sample IDs, chromosome names, and
expression-matrix column headers. Stray leading/trailing whitespace in these
fields is a common artefact of hand-edited tables and Excel exports, and it
would otherwise **silently break exact identifier matching** (e.g. a user's
`sample_columns = c("con1", "con2")` against matrix column names that carry a
trailing space). Each value is therefore normalized once on import; the package
performs no further string cleaning downstream.

This is a safe transformation: `trimws()` removes only leading/trailing
whitespace and never alters internal content, and it is applied solely to
identifiers (gene symbols, chromosome names, sample IDs), never to free text.
Genomic identifiers cannot legitimately encode surrounding spaces as part of
their identity, so normalization cannot change a valid value. It also only
drops entirely-blank entries, which are garbage input by any standard.

## 论证 / 实测依据

### 全库 63 处 trimws 分类

| 类别 | 数量 | 理由 | 判断 |
|---|---|---|---|
| 基因符号/分号拆分值 | ~55 处 | 用户文件基因名、`strsplit(x, ";")` 拆分后、SYMBOL 注解——外部数据，Excel 粘贴产生空格是真实场景 | 合理 |
| 列名/表头 | 3 处（analysis.R:270/1927、utils.R:1864）+ 表头 token 2 处 | 表达式矩阵列头来自外部工具/手编表格 | 合理 |
| 非空校验 | analysis.R:299/336、data_processing.R:496 | `nzchar(trimws(x))` 判断"全空白"字符串 | 合理 |

### 审稿人引用的那处已有注释（analysis.R:270-271）

```r
# Validate sample column names (trim before checking)
# Whitespace in headers is a common artefact of hand-edited tables/Excel exports.
sample_cols <- trimws(colnames(tpm_mat_raw))
```

直接回答审稿人的反问：**有理由**——`sample_columns` 要与 `colnames` 精确匹配，
带尾随空格会静默匹配失败/选错列。

### 无危害论证

- `trimws()` 只去首尾空白，不碰内部内容。
- 只作用于标识符（基因名/染色体名/样本 ID/表头），不作用于自由文本。
- 标识符在生物学上不可能以首尾空格作为身份一部分——规范化不会改变任何合法值。
- 唯一"副作用"是把纯空白输入滤掉，这是正确的垃圾输入处理。
- 输入边界一次性规范化，下游无额外字符串清洗。

## 战略定位

- **保留，零改动**：有正当理由、无危害、注释已就位。审稿人问的是"有没有理由"
  而非"有没有危害"——用一句话说明理由即可，不为回复而改代码。
- **`toupper()` 重复调用（审稿人 #10）已收敛**：全库逐函数排查后，修掉了 4 处
  **同一变量被重复 `toupper`** 的真实重复（`run_lfc_violin` 的 `valid_targets`、
  `run_go_enrichment` 的 `names(universe_genes)`、`run_heatmap_and_connectivity`
  的 `loop_stats_df[[gene_col_name]]` 与 `stats_subset$Gene`）——改为"提一次、
  子集/重排跟着走"（toupper 幂等，结果等价）。**补充：细化的两个闭包函数
  （`.refine_compute_targets` ×2、`.refine_target_annotations` ×3）里
  `toupper(whitelist)` 共 5 处重复调用也已收敛**——各函数入口加
  `whitelist_upper <- toupper(whitelist)`，闭包捕获复用（`whitelist` 是参数且
  未修改，安全）；另修 `.refine_load_validate_data` 的 `toupper(anno_genes)`
  ×2。其余 ~50 处是**不同变量的输入规范化**，无法也无需提。覆盖测试
  （test-analysis-modules 44 / test-coverage-analysis 46 / test-mock-modules 23
  / test-refine 29 / test-review-fixes 250）0 失败。
- **`log_message()` 工厂**：全局定义一次 + 8 函数取 quiet 闭包，设计合理，不改。

---

# Part 6: `.require_pkg()` 统一可选包守卫（审稿人"commonly repeated tasks"最对口项）

## 背景

审稿人总建议点名 "define accessory utils functions for commonly repeated tasks
and use these functions wherever possible"。全库 `if (!requireNamespace(...))`
守卫有三种语义各写各的、消息措辞并存（`"install 'TFBSTools' package."` vs
`"' not installed"` vs `"install 'networkD3' and 'htmlwidgets' packages."`）——
教科书级的重复任务。

## 修复

新增 `.require_pkg(pkg, feature = NULL, on_missing = c("stop", "warn", "return"))`
（`R/utils.R`），统一消息模板：

> `Package 'X' is required for <feature>. Install with BiocManager::install('X').`
> （warn/return 模式追加 " Skipping."）

**迁移 17 处**纯守卫（stop 9 / warn 2 / message+return 5），按语义映射：
- `on_missing = "stop"`：硬依赖（org_db、STRINGdb、ggraph、rmarkdown、txdb/orgdb、
  rtracklayer×2、可视化 arg/pkg）
- `on_missing = "warn"`：警告后由调用方分支（JASPAR2020 motif analysis/enrichment）
- `on_missing = "return"`：消息后由调用方返回各自值（TFBSTools/JASPAR2020/ggseqlogo）
- 返回 `invisible(TRUE)`（已装）/ `invisible(FALSE)`（缺失），调用方 `if (!.require_pkg(...))` 分支

**刻意保留的裸调用（4+1 处）**：
- 4 处组合条件（motifmatchr、networkD3‖htmlwidgets×2、ComplexHeatmap）——非单包
  守卫，helper 表达不了，保留裸调用（混合设计，符合 "wherever possible"）
- 1 处 `annotation.R:4639`（rtracklayer → **warn 后继续 fallback 到 BED-only**）——
  "warn+继续"语义匹配不了 stop/warn/return 任何一档，保留原样并加注释说明

## 稳定性验证

- **全量测试 1346 断言，0 失败，0 错误**（含 test-coverage-boost:1480 锁定的
  "is required for GO analysis"——统一模板仍含该子串，测试未改即通过）
- 三模式冒烟：stop 含 install 提示 / warn 返回 FALSE / return 返回 FALSE / 已装返回 TRUE
- 消息统一后所有调用点上下文经 `feature` 参数保留（如 "GO analysis"、
  "PPI analysis"、"Motif family annotation"）

## 战略定位

- **直接回应审稿人总建议**：这是 "commonly repeated tasks" 最典型的落地，从三种
  措辞并存 → 一个 helper + 统一模板。
- **边界诚实**：只覆盖 17/22 纯守卫；组合条件与"warn+继续"语义如实保留并说明——
  不假装"全覆盖"，审稿人看重的是"用到位"而非"全替换"。

---

# Part 5: README `eval=FALSE` 口径修正（第一轮遗留，内部审查发现）

## 背景

第一轮回复称 "Only two README chunks remain intentionally unevaluated:
installation and `looplook_report()`"。内部审查发现 README 实际有 **3 个**
`eval=FALSE`——多出的 L214（`validate_epeG_by_chromatin` /
`refine_loop_anchors_by_chromatin` 示例）用的是**不存在的占位路径**
（`"H3K4me1_peaks.bed"`，非 `system.file`），既不能运行也没写说明，
与第一轮口径矛盾。

## 修复（README.Rmd L214）

- 占位路径 → **真实打包示例**：`system.file("extdata", "example_h3k4me1_peaks.bed", ...)`
  等（extdata 本就有这些文件）。
- `eval=FALSE` → `eval = requireNamespace("looplook", quietly = TRUE)`，
  与其他模块 demo chunk 一致。
- 补注释：示例用打包染色质文件，用户需替换为自己的 ChIP-seq/CUT&Tag/ATAC-seq 文件。

## 实测（示例数据上串起 README 全流程）

`refined_res` 在 README 构建中已由 L180（`refine_loop_anchors_by_expression`，
`eval = requireNamespace`）计算——L214 引用的是**合法变量**。实测：

- `validate_epeG_by_chromatin`：~0.1 s
- `refine_loop_anchors_by_chromatin`：~1.4 s

README 构建开销可接受，且该 chunk 现在真实可运行。

## 结果

README 的 `eval=FALSE` 回到 **2 个**（L73 安装、L319 looplook_report），
与第一轮回复口径**完全一致**。回复中补一句说明即可：

> We also verified the README's unevaluated chunks: one Module 3B demo had slipped
> in with placeholder file paths. It now uses the packaged example chromatin BED
> files and evaluates like the other module demos, so the README again has only the
> two intentionally unevaluated chunks (installation and `looplook_report()`).

## 战略定位

- **把声称做成真**：第一轮"只剩 2 个"的断言现在属实（之前是口径漏洞）。
- 教训：任何"数量性"声称都要逐项核对，不能凭记忆——这次是内部审查抓出来的。

---

# Part 7: R CMD check 计时（审稿人第一轮"10 分钟内"要求，实测+修复）

## 实测结果（本机 macOS arm64，R 4.5.x）

`R CMD check --no-build-vignettes` 实测 **21.71 分钟** —— 超出 10 分钟线两倍多。
0 ERROR / 0 WARNING / 1 NOTE（"future file timestamps"＝沙箱时钟伪影，非真实）。

耗时分解：

| 环节 | 耗时 |
|---|---|
| testthat 套件（`testthat.R`） | 9–13 分钟（check 环境，慢于 `test_dir`） |
| vignette 代码执行 | 147s（`--no-build-vignettes` 不跳过此项） |
| examples | ~40s（`looplook_report` 8.4s、`annotate_peaks_and_loops` 5.4s） |
| 其余检查 | ~1 分钟 |

## 瓶颈定位（逐 test_that 计时）

`test-analysis-modules.R` 一条测试独吃 **994s（16.6 分钟）**：

```
run_ppi_analysis uses explicit ppi_species_id   994.1s
run_go_enrichment with NULL universe_genes      21.5s
run_go_enrichment returns list with result       20.7s
run_go_enrichment handles unmapped genes         12.5s
其余 27 条                                            合计 ~3s
```

## 根因

该测试是 STRINGdb **网络测试**，但 gate 只有 `skip_if(is.null(species_id))`
（网络可达性）——构建服务器有网络即必跑，STRINGdb 首用下载物种映射表 → 16.6 分钟。
而同文件另一条网络测试（"External STRINGdb integration test"）用的是
`LOOPLOOK_RUN_NETWORK_TESTS=true` 环境变量 opt-in。**两条网络测试 gate 方式不一致**，
且回复文档 Part 6"网络测试 opt-in"的声称对这条不成立。

## 修复（已做）

给 L159 测试加与另一条一致的 env gate：

```r
skip_if(tolower(Sys.getenv("LOOPLOOK_RUN_NETWORK_TESTS", unset = "false")) != "true",
  "Network test disabled; set LOOPLOOK_RUN_NETWORK_TESTS=true to run.")
```

实测：test-analysis-modules.R **1053s → 69.4s**（43 passed / 0 failed / 2 skipped，
两条网络测试均按 env 开关跳过）。

## 修复后实测（重跑）

- **21.71 → 11.76 分钟**（带 `--as-cran`，0 ERROR / 0 WARNING / 1 NOTE 时钟伪影）。
- **cran=FALSE（builder 近似口径）：11.35 分钟，Status OK，0 ERROR / 0 WARNING / 0 NOTE。**
- testthat.R **224s**（原 9–13 分钟）、vignette 代码 143s、examples 31s。
- `--as-cran` 只多占 ~0.4 分钟；builder 口径下仍**超线 ~1.35 分钟**。
- 注意：本机是 macOS，Bioconductor Linux builder 耗时可能不同；11.35 是保守代理值。

## 决策点（已执行方案 3 并验证 → 无效 → 回退 → 如实上报）

1. ~~如实上报约 12 分钟~~
2. ~~再压一次冲进 10 分钟~~（砍 vignette 重复块 / 减慢测试 → 丢演示与覆盖）
3. **方案 3 实测无效**：`annotation_integrated` 的 target_bed 减采样 300→100 peak 后，
   chunk 耗时 **36.7 → 36.6s（无变化）**。这 36s 是**固定开销**（TxDb 加载 +
   ChIPseeker 全基因模型注释 + 图谱/枢纽 + 绘图），与 peak 数量无关。之前的
   "省 40s"预期基于线性假设，是错的。**已回退该改动**（避免无收益噪音 +
   "0/300"→"0/100" 警告文字变化 + 一条不实注释）。

**最终结论：如实上报。** R CMD check `--no-build-vignettes` 实测 **~11.3 分钟**
（cran=FALSE，builder 近似口径，0 ERROR / 0 WARNING / 0 NOTE）。超 10 分钟线
~1.3 分钟。为本机 macOS 代理值；Bioconductor Linux builder 耗时可能不同。
不为此砍 vignette 演示内容或测试覆盖（审稿人会反过来问为什么内容变少）。

## 教训

网络/重计算测试必须显式 opt-in，不能以"网络可达"为 gate——构建服务器永远可达。
另：更新回复里早先"约 3 分钟"的过期声称（那是加测试之前的数字）。

