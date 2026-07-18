# looplook R24 深度算法与代码审查报告

> **审查对象：仅本次上传的 `R24.zip`。**  
> 审查视角：生物信息学、表观遗传学、3D 基因组学、R/Bioconductor 软件工程。  
> 本报告重点回答两个问题：
>
> 1. 算法的数据流、状态更新和下游 API 是否形成闭环；
> 2. 生物学归因和代码实现是否足够精准，是否仍存在会改变结果的边界问题。

---

# 一、审查范围与限制

## 1. R24 源码结构

R24 包含：

```text
R/
├── analysis.R
├── annotation.R
├── data_processing.R
├── looplook-package.R
├── utils.R
└── visualization.R
```

源码规模：

| 文件 | 行数 |
|---|---:|
| `annotation.R` | 7,426 |
| `analysis.R` | 2,467 |
| `utils.R` | 2,085 |
| `data_processing.R` | 1,540 |
| `visualization.R` | 995 |
| `looplook-package.R` | 71 |
| **合计** | **14,584** |

R24 的六个 R 源文件与上一版 R23 的对应源文件一致；ZIP 哈希差异来自压缩包本身，不代表 R 代码发生变化。

## 2. 本轮执行的静态检查

已检查：

- 函数结构和调用关系；
- 字符串、注释和括号配对；
- loop consolidation；
- canonical anchor 构建；
- ChIPseeker/TxDb 注释；
- promoter/gene-body/distal 分类；
- topology 与多跳 target assignment；
- expression refinement；
- chromatin refinement；
- TSS 重注释；
- target/profile API；
- strict/filled fallback；
- long-table provenance；
- downstream LFC/GSEA/GO/heatmap/motif/PPI 逻辑；
- 文档与实现的一致性。

未发现明显的括号、字符串或函数截断问题。

## 3. 尚不能完成的验证

当前环境没有 R/Rscript，且 ZIP 中没有：

```text
DESCRIPTION
NAMESPACE
tests/testthat/
vignettes/
man/
.github/workflows/
```

因此本报告不能替代：

```text
R CMD check
BiocCheck
testthat
真实TxDb/OrgDb运行
真实ChIPseeker输出验证
真实bigWig导入
真实HiChIP数据回归测试
```

所有“代码正确”结论均限定为静态源码层面。

---

# 二、总体结论

## 1. 是否闭环

### 数据流与软件 API：基本闭环

当前主流程已经形成稳定闭环：

```text
BEDPE
→ canonical anchors / loop graph
→ positional P/G/E annotation
→ structural all/promoter target assignment
→ expression current-target filtering
→ chromatin-aware anchor/TSS refinement
→ strict/filled current target reconstruction
→ target_gene_links provenance
→ profile_target_genes()
```

尤其以下接口已经稳定：

```text
loops/all
→ Putative_Target_Genes

loops/promoter
→ Promoter_Target_Genes

targets/all/strict
→ Assigned_Target_Genes

targets/all/filled
→ Assigned_Target_Genes_Filled

targets/promoter/strict
→ Regulated_promoter_genes

targets/promoter/filled
→ Regulated_promoter_genes_Filled
```

Expression 是否已过滤，由传入对象的 refinement 阶段决定，不再通过 `active` mode 再切换。

### 生物学推断：尚未完全闭环

当前仍有几处会影响基因归属或 chromatin identity 的高优先级假设：

1. basic annotation 中未测量基因被当成表达量 0；
2. absolute expression threshold = 0 时，零表达基因会被当作 active；
3. G/eG→P 或 dual 后未强制验证 TSS；
4. H3K4me1/H3K4me3 bigWig 比值使用不同覆盖区域的 covered mean，且 NA ratio 默认走 promoter 分支；
5. 多跳 target assignment 中 promoter 身份可压过更近的 gene-body link。

因此：

> **R24 的工程与状态更新主链已经闭环，但要达到方法学论文或高标准 Bioconductor 方法包的“科学推断闭环”，仍需收紧上述边界。**

## 2. 严重度汇总

| 等级 | 数量 | 判断 |
|---|---:|---|
| P0 运行阻断 | 0 | 静态检查未发现 |
| P1 高优先级精准性问题 | 4 | 建议冻结前修复 |
| P2 方法学/边界问题 | 7 | 建议文档化或增加参数 |
| P3 工程与文档问题 | 若干 | 不直接改变核心结果 |

---

# 三、当前已经做得较好的部分

# 1. Loop consolidation 逻辑较成熟

`data_processing.R` 已实现：

- BEDPE 坐标和 anchor 顺序规范化；
- 0-based BEDPE → 1-based GRanges；
- score 自动识别及 p-value-like 警告；
- consensus/intersect/union 三种模式；
- replicate-balanced score 聚合；
- connected-component chaining 诊断；
- blacklist/ROI post-filter；
- provenance metadata。

特别是 consolidation 阶段已经明确说明：

```text
intersect mode不是对称操作
connected-component可能产生transitive chaining
warn模式不会自动删除宽cluster
```

这比许多常见 loop merge 脚本更透明。

# 2. Positional annotation 与 functional identity 已分开

R24 明确将：

```text
P = promoter-positioned
G = gene-body-positioned
E = distal/intergenic-positioned
```

而不是直接把 E 等价为功能 enhancer。

这是正确的表观遗传学边界。

# 3. Expression refinement 的状态模型已经合理

当前能够区分：

```text
active
measured_silent
unmeasured
measured_not_assessed
```

并保留：

```text
eP:GeneA
eG:GeneB
```

而不是删除 positional gene。

Loop summary 使用 current target set，long-table 保存完整 structural links，设计正确。

# 4. Chromatin mark 使用三态语义

每个 mark 区分：

```text
TRUE  = peak detected
FALSE = BED提供但此处未call peak
NA    = mark未提供
```

代码和文档也承认：

```text
FALSE ≠ 生物学绝对缺失
```

这点十分重要。

# 5. Chromatin conflicting state 比较保守

H3K27me3 与 active marks 同时存在时，代码输出：

```text
conflicting_marks
```

并保留原 positional type，不直接断言为高置信度 active promoter、enhancer 或 dual。

这避免了把：

- 真正 bivalent；
- 细胞群体异质性；
- peak boundary overlap；
- broad-domain overlap；

混为一个确定生物学状态。

# 6. Target/profile 当前主链已经统一

Chromatin finalization 现在采用：

```text
all structural links
→ expression-qualified current_links
→ current promoter-first priority
→ strict target columns
→ current linear fallback
→ current promoter evidence
→ membership flags
→ retained_after_refinement
```

而完整 links 继续保留为 provenance。

这一部分已经闭环。

# 7. TSS alignment 的防错较强

`.reannotate_tss_genes()` 不再盲目信任行顺序，而使用：

```text
ID match
→ coordinate join
→ duplicate-coordinate ambiguity
→ unmatched → NA
```

匹配不到时宁可清空 gene，也不错误地按行对齐，是正确的保守策略。

---

# 四、P1 高优先级问题

# P1-1：Basic annotation 把“未测量”基因当成 TPM = 0

## 代码位置

`utils.R:791–801`

```r
gene_map$tpm <- if (!is.null(gene_expr_map)) {
    expr_upper <- setNames(gene_expr_map, toupper(names(gene_expr_map)))
    val <- expr_upper[toupper(gene_map$SYMBOL)]
    ifelse(is.na(val), 0, val)
} else {
    0
}
```

随后：

```r
gene_map$is_active <- ...
```

并在 `utils.R:666–675` 中：

```r
has_active = any(is_active)
filter(!has_active | is_active)
```

## 实际含义

当 expression matrix 不包含某个候选 gene 时：

```text
unmeasured
```

被直接变成：

```text
TPM = 0
measured_silent
```

只要同一 promoter overlap 中存在另一个 measured-active gene，未测量候选就会被删除。

## 示例

一个 anchor 同时靠近：

```text
GeneA：TPM = 5
GeneB：expression matrix中不存在
```

代码可能只保留：

```text
GeneA
```

但 GeneB 的真实状态是：

```text
unknown/unmeasured
```

而不是 silent。

这与后续 expression refinement 精心区分：

```text
measured_silent
unmeasured
```

的设计不一致。

## 生物学后果

该问题可能在以下场景中改变 basic gene assignment：

- expression matrix 只包含 protein-coding genes；
-低表达 genes 在预处理时被删除；
- lncRNA/pseudogene 不在矩阵中；
- SYMBOL/ENSEMBL 映射不完整；
-自定义基因集或非模式物种；
-单细胞 pseudo-bulk 基因过滤。

## 推荐修复

保留显式 measured 状态：

```r
val <- expr_upper[toupper(gene_map$SYMBOL)]

gene_map$measured <- !is.na(val)
gene_map$tpm <- as.numeric(val)

gene_map$is_active <- gene_map$measured &
    if (min_expr == 0) gene_map$tpm > 0 else gene_map$tpm >= min_expr
```

Conflict resolution 应提供明确策略：

```text
unmeasured_policy = "keep"      # 推荐默认
unmeasured_policy = "penalize"
unmeasured_policy = "drop"
```

推荐规则：

1. biotype/位置优先级照常；
2. measured-active 可以优先于 measured-silent；
3. 不应仅因缺少测量而删除 unmeasured candidate；
4. 如果 expression mapping coverage 过低，应禁用 expression-assisted conflict resolution。

## 严重度

```text
P1
```

它直接改变 basic annotation gene assignment。

---

# P1-2：Absolute threshold = 0 会把零表达基因判定为 active

## 代码位置

`annotation.R:2636`

```r
whitelist <- names(vals)[
    vals >= effective_threshold &
    !is.na(vals)
]
```

Link state 也使用：

`annotation.R:2169`

```r
passes <- Mean_Expression >= threshold
```

但 basic conflict resolution 在 `utils.R:798–801` 中使用：

```r
if (min_expr == 0) tpm > 0
```

## 问题

文档将：

```text
threshold = 0
```

描述为：

> retain any detected expression

但当前 refinement 中：

```text
expression = 0
```

也满足：

```r
0 >= 0
```

会被判定为 active。

因此 basic annotation 与 refinement 对 threshold 0 的定义不一致：

| 阶段 | threshold 0 |
|---|---|
| basic conflict | `> 0` |
| expression refinement | `>= 0` |

## 后果

用户希望保留“任何非零表达”时，零表达 gene 不会被过滤：

```text
Putative_Target_Genes
Promoter_Target_Genes
Assigned_Target_Genes[_Filled]
```

均可能保留零表达基因。

## 推荐修复

建立统一 helper：

```r
.passes_expression_threshold <- function(x, threshold, mode) {
    if (mode == "absolute" && threshold == 0) {
        return(!is.na(x) & x > 0)
    }
    !is.na(x) & x >= threshold
}
```

以下位置统一调用：

- whitelist 构建；
- target_gene_links；
- promoter stats；
- chromatin expression recovery；
- Has_Active_Target。

同时增加验证：

```text
absolute threshold必须有限且>=0
quantile threshold必须在[0,1]
至少存在一个有限expression value
```

## 严重度

```text
P1
```

默认 threshold = 1 不受影响，但 threshold = 0 是明确支持的公开用法。

---

# P1-3：G/eG→P 或 dual 后没有强制验证 TSS

## 相关代码

Chromatin reclassification：

`annotation.R:4494–4512`

```text
eG + promoter_like → P
G + H3K4me3 → P
G + H3K4me1/H3K4me3 high ratio → dual
```

Gene restoration：

`annotation.R:4815–4836`

当前只有：

```r
requires_tss <-
    old_type == "E" &
    new_type %in% c("P", "dual")
```

因此：

```text
E→P/dual
```

必须重新用 TxDb 检查 TSS。

但：

```text
G→P
eG→P
G→dual
eG→dual
```

会直接沿用原 gene-body gene，除非 gene 本身缺失。

## 为什么不够严谨

Intragenic H3K4me3 可能代表：

- 同一基因的 alternative promoter；
- 相邻或反义基因的 promoter；
-未注释转录本 TSS；
- broad H3K4me3 domain；
-转录相关信号扩展；
-混合细胞群体中的 peak overlap。

“处于 GeneA gene body”并不能自动证明该 H3K4me3 peak 是：

```text
GeneA promoter
```

代码注释将其解释为 alternate TSS，但没有实际验证 TSS。

## 可能后果

一个 G anchor 被改为：

```text
P:GeneA
```

随后会进入：

```text
Promoter_Target_Genes
Regulated_promoter_genes
promoter_centric_stats
```

即使该 anchor 附近没有 GeneA 的已注释 TSS。

## 推荐修复

将 mandatory TSS check 扩展到：

```r
requires_tss <-
    old_type %in% c("E", "G", "eG") &
    new_type %in% c("P", "dual")
```

更精确的处理：

### 有 TSS 支持

```text
final P/dual
gene = TSS-supported gene
```

### 没有 TSS 支持

不应直接将 host gene 作为 promoter gene。

可选择：

```text
保留原 G/eG type
chromatin_state = promoter_like_no_annotated_TSS
```

或：

```text
final type = P
gene = NA
TSS_supported = FALSE
```

前者在现有 schema 中更保守。

## 严重度

```text
P1
```

它可改变 promoter target 与 promoter connectivity 统计。

---

# P1-4：H3K4me1/H3K4me3 bigWig ratio 的数学定义和 NA fallback 不够稳健

## 代码位置

`annotation.R:6942–6986`

当前分别计算：

```text
H3K4me1 covered_mean
= H3K4me1 AUC / H3K4me1 covered width

H3K4me3 covered_mean
= H3K4me3 AUC / H3K4me3 covered width
```

随后：

```r
log2_ratio =
    log2(H3K4me1_covered_mean + pseudocount) -
    log2(H3K4me3_covered_mean + pseudocount)
```

## 问题一：两个分母不是同一 genomic support

H3K4me1 和 H3K4me3 的 covered width 可以不同。

例如：

```text
H3K4me1：窄而高
H3K4me3：宽而中等
```

covered-mean ratio 可能显示 H3K4me1 dominant，即使 H3K4me3 在整个 anchor 的总信号更大。

如果要比较同一 anchor 上两个 marks 的相对 dominance，更自然的是使用相同 denominator：

```text
anchor-wide mean
= AUC / anchor width
```

两个 mark 的 anchor width 相同，因此其 ratio 等价于：

```text
AUC_H3K4me1 / AUC_H3K4me3
```

当前注释称 anchor-wide denominator 会引入 anchor size confounding，但在同一 anchor 内做 ratio 时，anchor width 会抵消，这一数学解释不准确。

## 问题二：NA ratio 被默认为 promoter-dominant

`annotation.R:6987–6993` 明确提示：

```text
NA ratio → classified as promoter-dominant
```

而 `.chromatin_reclassify()` 中：

```r
is_true_dual = FALSE
```

会使 dual-positive anchor进入 P/G 分支。

但 NA 可能来自：

- bigWig 缺失该 chromosome；
- seqlevel 不匹配；
- signal track存在空洞；
-一个 mark 没有导入到 anchor；
- track normalization/coverage 问题。

缺少定量数据不等于 promoter dominance。

## 问题三：固定 threshold = 3 缺少样本内校准

H3K4me1 和 H3K4me3 bigWig 必须具有可比尺度：

- 同一 normalization；
-相近测序深度；
-相同 signal definition；
-无不同 global scaling。

否则固定 ratio threshold 无普适意义。

## 问题四：“true dual”表述过强

双 peak + ratio threshold 只能支持：

```text
dual-marked / H3K4me1-dominant chromatin signature
```

不能直接证明：

```text
dual promoter-enhancer function
```

## 推荐修复

### 定量指标

首选：

```r
mean1_anchor <- AUC1 / anchor_width
mean3_anchor <- AUC3 / anchor_width

log2_ratio <- log2(mean1_anchor + pc) -
              log2(mean3_anchor + pc)
```

或固定宽度窗口、同一 bin grid 上计算 signal。

### Threshold

使用 reference distributions 校准：

```text
known promoters
known enhancers
```

例如根据 ROC、mixture model 或 empirical quantiles 选择 cutoff，而不是固定 3。

### NA

```text
NA ratio → dual_like_unresolved / uncertain
```

不能默认为 promoter。

### 命名

将：

```text
is_true_dual
```

改为：

```text
is_me1_dominant_dual_signature
```

## 严重度

```text
P1
```

它直接影响 P/E/dual reclassification，是当前 chromatin 方法中最需要验证的部分。

---

# 五、P2 方法学与边界问题

# P2-1：多跳 target assignment 先看 promoter role，未优先考虑 path length

## 代码位置

Basic topology：

`annotation.R:1287–1307`

```text
loop ego order = neighbor_hop
target ego order = neighbor_hop + 1
```

Promoter priority：

`annotation.R:144–151`

```r
只要ego内存在promoter anchor
→ 只使用promoter genes
```

Chromatin current assignment：

`annotation.R:5289–5303`

```r
promoter priority = 1
gene_body priority = 2
```

没有把：

```text
path_length
```

纳入 priority。

## 风险场景

当 `neighbor_hop > 0`：

```text
direct 1-hop active gene-body gene
2-hop promoter gene
```

当前可能选择：

```text
2-hop promoter gene
```

而丢弃更近的 gene-body-supported gene。

这不是程序错误，而是隐含的方法选择：

```text
role-first
```

而不是：

```text
distance-first
```

## 推荐修复

建议默认采用 lexicographic priority：

```text
1. 最小 path_length
2. 同 path_length 内 promoter > gene_body
3. 同优先级内保留并列 genes
```

或增加：

```r
target_priority = c(
    "distance_then_role",
    "role_then_distance"
)
```

默认建议：

```text
distance_then_role
```

同时将：

```text
direct targets
expanded targets
```

在 summary 或 profile 中分层，不能只在 long table 保留 provenance。

## 严重度

```text
P2
```

默认 `neighbor_hop = 0` 时风险较低；扩展 hop 时影响明显。

---

# P2-2：`anchor_merge_gap > 0` 存在 transitive anchor chaining

## 代码位置

`annotation.R:1081–1086`

```r
GenomicRanges::reduce(
    gr_raw,
    min.gapwidth = anchor_merge_gap + 1
)
```

`reduce()` 是传递性合并：

```text
A接近B
B接近C
→ A/B/C合成一个anchor
```

即使 A 与 C 之间超过 `anchor_merge_gap`。

## 后果

可能形成较宽 canonical anchor，进而改变：

- P/G/E 分类；
- target BED overlap；
- chromatin mark overlap；
- loop self-collapse；
- graph degree；
- promoter connectivity。

Loop consolidation 阶段已有 chaining diagnosis，但 annotation anchor merge 没有对应 span guard。

## 推荐修复

增加：

```text
merged_anchor_width
n_raw_anchors_per_merged_anchor
max_pairwise_gap
span/raw_median_width ratio
```

并增加：

```r
anchor_merge_policy = c("warn", "drop", "error", "none")
max_merged_anchor_width
```

默认 `anchor_merge_gap = 0` 是安全的，应继续保持。

---

# P2-3：未知 ChIPseeker annotation 类型默认映射为 E

## 代码位置

`annotation.R:68–99`

文档声称可能返回：

```text
Unknown
```

但代码最终：

```r
"E"
```

## 风险

如果 ChIPseeker 新版本、custom TxDb 或特殊 feature label产生未知类别，代码会将其当作 distal/intergenic，进入 E-based loop type 和 target逻辑。

“未知”不等于“distal”。

## 推荐修复

```r
return("Unknown")
```

并规定：

```text
Unknown不参与P/G/E target assignment
```

同时在输出和 warning 中要求用户检查 mapping。

---

# P2-4：`anchor_gap` 在 mandatory physical-overlap 规则下并不是真正的 proximity tolerance

## 代码位置

Target overlap：

`annotation.R:226–250`

Chromatin mark overlap：

`annotation.R:7281–7330`

流程是：

```text
findOverlaps(maxgap = anchor_gap)
→ 再要求 physical overlap >= anchor_min_overlap
```

而：

```text
anchor_min_overlap >= 1
```

是强制验证条件。

因此不重叠但距离 anchor 200 bp 的 peak：

```text
即使anchor_gap = 200
```

也会因为 physical overlap = 0 被删除。

## 判断

当前 `anchor_gap` 更像候选 hit 预筛选参数，而不是最终 gap tolerance。

文档中“gap tolerance”容易让用户误以为 proximity-only peak 会被接受。

## 推荐方案

二选一：

### 严格 overlap 模式

删除或弱化 `anchor_gap`，明确只接受真实 overlap。

### Proximity 模式

允许：

```text
anchor_min_overlap = 0
```

并新增：

```r
overlap_mode = c("physical", "proximity")
```

避免一个参数同时表达两种含义。

---

# P2-5：Hub 主要按 raw loop rows 判定，易受重复/平行边影响

## 代码位置

Promoter stats：

`annotation.R:1419–1479`

Distal stats：

`annotation.R:1541–1588`

虽然同时输出：

```text
Total_Loops
n_Unique_Contacts
```

但 high-connectivity cutoff 使用：

```text
Total_Loops
```

`Total_Loops` 包含 anchor merge 后的 parallel loop rows。

## 后果

如果输入中：

-重复 loop rows；
-多个 callers 的重复输出；
-样本支持以重复行形式保留；
-未先 consolidation；

某些 promoter 或 distal anchor 的 hub status 会被放大。

## 推荐修复

增加：

```r
hub_metric = c(
    "unique_contacts",
    "total_loops",
    "support_weighted"
)
```

默认建议：

```text
unique_contacts
```

Replicate support 单独作为 confidence，不应和 topology degree 混为一个指标。

---

# P2-6：Expression refinement 允许在 refined object 上再次运行

## 代码位置

`annotation.R:3449–3462`

当前只 warning，随后继续执行。

## 问题

第一次 refinement 已经：

-过滤 current target columns；
-修改 P/G 为 eP/eG；
-写入 expression state。

再次 refinement，尤其使用更低 threshold 时，无法恢复第一次已删除的 current genes。

因此 sensitivity analysis 结果取决于调用历史，而不是单一 threshold。

## 推荐修复

默认 stop：

```r
allow_rerefine = FALSE
```

若确实允许：

```text
只能从保存的 original/basic state 重建
```

不能对 current filtered object 再过滤。

---

# P2-7：Profile 输入和统计背景需要更强验证

## 代码位置

`analysis.R:217–230`

Profile 直接：

```r
read_robust_general()
rowMeans(tpm_mat_raw, na.rm = TRUE)
```

没有复用 `load_expression_matrix()` 中的：

- non-numeric token validation；
- duplicate gene处理；
- duplicate sample检查；
- case-collision诊断；
-全 NA 行处理。

## 具体风险

### Expression matrix

字符 token、重复 gene ID 或 case collision 可能造成：

- `rowMeans()` 报错；
-重复 gene被隐式处理；
-heatmap exact matching失败；
- LFC、GO 和 heatmap使用不同匹配逻辑。

### GO universe

GO background 使用：

```text
diff_file中具有finite LFC的全部genes
```

只有在 `diff_file` 是完整 DESeq2/edgeR tested-gene table 时才合理。

如果用户只提供 DEG 子集：

```text
GO universe严重偏倚
```

### Heatmap

部分位置使用 exact symbol matching，而 LFC 模块使用 case-insensitive matching，可能导致同一 gene 在一个模块被识别、另一个模块丢失。

## 推荐修复

1. Profile 复用统一 expression parser；
2. 明确要求 `diff_file` 包含全部 tested genes，而非仅 DEG；
3.允许独立提供：

```r
universe_genes
```

4.统一 identifier normalization/mapping；
5.检查 metadata SampleID 唯一且 Group 非空。

---

# 六、其他需要明确的生物学解释边界

# 1. Physical loop 不等于调控因果

以下列：

```text
Putative_Target_Genes
Promoter_Target_Genes
Assigned_Target_Genes
Regulated_promoter_genes
```

均代表：

```text
loop/topology-supported candidate relationship
```

不能直接证明：

```text
该element调控该gene
```

建议在文档中反复使用：

```text
candidate
loop-supported
putative
```

而不是 causal wording。

`Regulated_promoter_genes` 名称本身偏强。为保持 API 不破坏，可以保留列名，但必须明确其 operational definition。

# 2. H3K4me3 不等于已验证 promoter function

H3K4me3 是强 promoter-associated mark，但：

- broad domain；
- cryptic transcription；
- alternative TSS；
- eRNA TSS；
-混合细胞群体；

均可能产生 H3K4me3 signal。

因此 H3K4me3+ 应解释为：

```text
promoter-like chromatin
```

而非单独证明 promoter-gene relationship。

# 3. H3K4me1/H3K4me3 双阳性不证明 dual function

双阳性及 ratio 只能定义 chromatin signature。

功能上的 promoter-enhancer duality 仍需要：

- CAGE/PRO-cap/NET-CAGE；
- reporter/CRISPRi；
- perturbation RNA-seq；
- allele-specific/eQTL；
- nascent transcription；
- enhancer RNA/TSS directionality；

验证。

# 4. `conflicting_marks` 不应直接称为 bivalent

当前代码已经较保守。

建议结果表述为：

```text
conflicting/bivalent-like chromatin signature
```

而非确定的 single-cell bivalency。

# 5. Expanded topology targets应作为较低层级证据

建议 evidence hierarchy：

```text
Tier 1: local/direct promoter link
Tier 2: direct gene-body link
Tier 3: one-hop expanded promoter
Tier 4: multi-hop expanded target
Tier 5: linear fallback
```

不要在主结果中无分层合并。

---

# 七、各阶段闭环评价

| 阶段 | 数据流闭环 | 生物学精准性 | 评价 |
|---|---:|---:|---|
| BEDPE import | 高 | 高 | 坐标与score处理较成熟 |
| Loop consolidation | 高 | 中高 | chaining已诊断，但仍依赖参数 |
| Canonical anchor merge | 中高 | 中 | >0 gap有transitive merge风险 |
| Basic P/G/E annotation | 高 | 中高 | unmeasured=0问题需修复 |
| Basic target assignment | 高 | 中高 | hop expansion和priority需分层 |
| Expression refinement | 高 | 高 | threshold=0语义需统一 |
| Chromatin validation | 高 | 中高 | mark三态优秀；ratio需重审 |
| TSS/gene restoration | 高 | 中 | G/eG→P缺mandatory TSS |
| Target strict/filled | 高 | 高 | R23/R24已闭环 |
| Long-table provenance | 高 | 高 | current/full关系清晰 |
| Profile API | 高 | 中高 | 输入与universe需强化 |
| Causal interpretation | 不适用 | 中 | 必须保持candidate定位 |

---

# 八、建议的最小修复顺序

# 第一优先级：修复 expression 语义

1. 未测量 gene 不再自动赋 0；
2.统一 `threshold = 0` 为 `expression > 0`；
3.增加 absolute/quantile threshold validation；
4.建立单一 `.passes_expression_threshold()` helper。

# 第二优先级：收紧 promoter gene attribution

1. G/eG→P/dual 强制 TSS check；
2.无 TSS 支持时不要保留 host gene作为 promoter gene；
3.输出明确的 `TSS_supported` / `promoter_like_no_TSS` 状态。

# 第三优先级：重新验证 bigWig ratio

1.改为同一 genomic denominator；
2. NA ratio → unresolved；
3. threshold由 reference promoter/enhancer校准；
4.降低“true dual”措辞。

# 第四优先级：收紧 topology priority

1. path length优先；
2. promoter role作为第二层 priority；
3. direct/expanded target分层；
4. neighbor_hop > 0 时输出 expansion diagnostics。

# 第五优先级：工程与统计健壮性

1. anchor merge chaining diagnostics；
2. hub metric可配置；
3. profile统一输入解析；
4.独立 GO universe；
5.禁止在 refined object 上直接重复 refine。

---

# 九、正式测试矩阵

## A. Basic annotation

### Test 1：active + unmeasured candidate

```text
GeneA active
GeneB unmeasured
```

不得自动把 GeneB 当 silent 删除。

### Test 2：all candidates unmeasured

应保留不确定候选或明确标记 unresolved。

### Test 3：unknown annotation category

应为：

```text
Unknown
```

不能自动变 E。

---

## B. Expression

### Test 4：absolute threshold = 0

```text
expression 0 → measured_silent
expression >0 → active
```

### Test 5：quantile threshold边界

```text
threshold <0或>1 → stop
```

### Test 6：all NA expression

明确 stop 或返回 not_assessed，不得静默产生空 whitelist。

### Test 7：repeat refinement

默认 stop。

---

## C. Chromatin/TSS

### Test 8：G→P，有同gene TSS

保留 TSS-supported gene。

### Test 9：G→P，无任何TSS

不得直接输出：

```text
P:host_gene
```

### Test 10：G→P，附近TSS属于另一个gene

必须使用 TSS-supported gene。

### Test 11：E→P，无TSS

gene保持 NA，且 target/provenance一致。

---

## D. bigWig

### Test 12：相同AUC、不同covered width

验证 covered mean ratio和anchor-wide ratio差异。

### Test 13：某mark bigWig缺失

应输出 unresolved，不得自动 promoter。

### Test 14：不同global scaling

要求校准或警告。

---

## E. Topology

### Test 15：1-hop gene-body + 2-hop promoter

验证 distance-first policy。

### Test 16：同hop promoter + gene-body

promoter可优先。

### Test 17：neighbor_hop=0

保持当前保守默认结果。

---

## F. Profile

### Test 18：DEG-only diff file

应 warning universe不完整。

### Test 19：expression字符值/重复gene/sample

明确报错。

### Test 20：case mismatch

所有模块应使用一致的 ID mapping。

---

# 十、评分

## 1. Target schema / Profile API

```text
9.1 / 10
```

该部分已经闭环。

## 2. Basic annotation

```text
8.1–8.4 / 10
```

主要扣分来自 unmeasured=0 和 gene conflict语义。

## 3. Expression refinement

```text
8.6–8.8 / 10
```

状态设计优秀；threshold=0与输入验证需修复。

## 4. Chromatin refinement

```text
8.0–8.4 / 10
```

三态 mark和provenance较强；TSS与bigWig ratio仍是关键问题。

## 5. Topology target assignment

```text
8.2–8.5 / 10
```

默认 hop=0较稳健，多跳 priority需收紧。

## 6. Downstream profiling

```text
7.8–8.2 / 10
```

适合作为 exploratory pipeline，但输入、universe和统计匹配需加强。

## 7. 当前源码整体

```text
8.3–8.6 / 10
```

已经是高质量研究代码，但尚不能称为方法学层面完全冻结。

## 8. Bioconductor-ready

```text
5.5–6.0 / 10
```

主要限制：

```text
无DESCRIPTION/NAMESPACE
无testthat
无vignette
无CI
无R CMD check/BiocCheck
无真实依赖环境回归
```

---

# 十一、最终判断

## 是否存在明显代码逻辑崩坏

没有。

R24 未发现新的 P0 或 target/profile 主链回归。

## Target/profile API 是否闭环

是。

以下维度已经稳定：

```text
all vs promoter
strict vs filled
basic vs expression-refined vs chromatin-refined object
current summary vs full provenance links
```

## 整个算法是否可以宣称完全闭环

还不能。

当前剩余的核心问题不是字段或 API，而是更深层的生物学归因：

```text
unmeasured gene如何处理
zero expression如何定义
intragenic H3K4me3如何绑定TSS/gene
H3K4me1/H3K4me3 ratio如何定量
multi-hop target如何排序
```

这些问题修复后，才能更有把握地称为：

```text
算法定义闭环
实现闭环
生物学语义闭环
```

## 最终一句话

> **R24 的软件架构、refinement 状态机、target/profile API 与 provenance 已经闭环；但 basic expression conflict、G/eG→P 的 TSS 验证和 dual-mark bigWig ratio 仍会实质影响生物学归因。建议暂时冻结现有 API，不再改字段结构，只针对这三个算法核心做最后一轮修订与正式测试。**
