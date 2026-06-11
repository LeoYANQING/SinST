# 空间组学数据处理标准流程

> 给自己做的速查笔记。第一次拿到 ST `.h5ad` 数据时，按这个顺序往下走。
> 关联：项目里的 ARTISTA 数据已经做完了 ①~⑨ 步，所以日常使用只需要 inspect + 可视化。

---

## TL;DR

完整流程是 **8 步常规 + 3 步 ST 特有**，**不是** "inspect → PCA → UMAP" 三步。

```
原始 .h5ad
   │
   ▼
 ① Inspect       看 shape / obs / var / X / layers
 ② QC            过滤低质量细胞和基因
 ③ Normalize     每个细胞 counts → 同一 scale
 ④ Log1p         log(x+1) 把长尾分布压扁
 ⑤ HVG           选 top 2000 高变基因
 ⑥ PCA           27000-d → 50-d
 ⑦ Neighbors     PCA 空间建 KNN 图
 ⑧ UMAP/Leiden   2-d 可视化 + 聚类
 ─── 以下 ST 特有 ───
 ⑨ Spatial neighbors      基于 obsm['spatial']，不是 ⑦ 的表达 KNN
 ⑩ Spatial autocorrelation  Moran's I：找空间上有 pattern 的基因
 ⑪ Niche / 区域分析          cell type co-localization, MENDER 微环境
```

关键洞察：PCA 和 UMAP 是工具，**真正决定结果的是 ②~⑤**（QC + 标准化 + HVG）。新手最常犯的错就是跳过这几步直接 PCA。

UMI 是什么：随机标签
测序前，每条 mRNA 被贴上一个随机 8~12 碱基的标签（Unique Molecular Identifier）：
  原始 mRNA:    GENE_A
                GENE_A
                GENE_A
                GENE_B

  贴 UMI 后:    GENE_A + UMI=AGGT
                GENE_A + UMI=CTGA
                GENE_A + UMI=TGCA
                GENE_B + UMI=GGTC

  扩增 + 测序 (每条复制 1000×):

     读到 1000 条:  GENE_A + AGGT     ┐
     读到 1000 条:  GENE_A + CTGA     │  原始 3 条 GENE_A
     读到 1000 条:  GENE_A + TGCA     ┘
     读到 1000 条:  GENE_B + GGTC     ──  原始 1 条 GENE_B

  去重计数:
     GENE_A: 看到 3 个不同的 UMI  →  原始有 3 条 mRNA
     GENE_B: 看到 1 个不同的 UMI  →  原始有 1 条 mRNA

  UMI counts = 不同的 UMI 标签数量 ≈ 原始 mRNA 分子的真实数量

  (没有 UMI 之前用 read counts，会被 PCR 扩增偏差严重污染。)

total UMI = 3 + 1 + 0 + ... + 2 = 这个细胞总共测到了多少条 mRNA
总 UMI ≈ 这个细胞"被测得有多好"
正常都是1000-30000
大于50000就是两个细胞黏在一起

adata.layers['counts'] = adata.X.copy()

---

## 每步在解决什么问题

| 步骤 | 解决什么问题 | 不做会怎样 |
|---|---|---|
| ① Inspect | 知道你有什么数据 | 后面瞎用 |
| ② QC | 过滤掉死细胞 (mt% 高)、空 droplet (基因数太少)、双胞 (基因数太多) | PC1 = "细胞是死是活" |
| ③ Normalize | 不同细胞测到的 UMI 总数差 10×，要拉到同一 scale | PC1 = "测序深度" |
| ④ Log1p | 表达量是长尾分布 (少数基因极高)，log 后接近高斯 | 距离被高表达基因主导 |
| ⑤ HVG | 27324 个基因里大部分是恒定噪声 | 信号被稀释，UMAP 一团糊 |
| ⑥ PCA | 线性降维，去冗余，加速下游 | UMAP / KNN 在 27000-d 上算极慢且不稳 |
| ⑦ Neighbors | 用 KNN 图描述细胞间相似性 | UMAP / 聚类没法跑 |
| ⑧ UMAP | 2-d 可视化；Leiden 在图上聚类 | 你看不见 |
| ⑨ Spatial neighbors | 描述细胞间**物理**相邻关系 | 没法做空间分析 / GNN |
| ⑩ Moran's I | 量化基因表达的空间结构性 | 找不到 spatially-variable genes |
| ⑪ Niche | 量化细胞类型的空间共现 | 错过组织区域结构 |

---

## ST 比 scRNA-seq 多 3 步：为什么？

scRNA-seq 数据只关心"细胞像谁"——所有距离都在表达空间里算。

ST 还关心"细胞挨着谁"——必须同时用 `obsm['spatial']`（物理坐标）和 `adata.X`（表达）两套距离。

⚠️ **不要混淆**：
- `obsp['connectivities']` = 表达空间 KNN (第 ⑦ 步产物)
- `obsp['spatial_connectivities']` = 物理空间 KNN (第 ⑨ 步产物)

stGraphFM / stVCR / GraphST 这类方法用的是**后者**。

---

## ARTISTA 已经处理到哪一步

BGI 团队 release `.h5ad` 时已经跑完前 ⑨ 步。验证字段对应：

| 字段 | 流程位置 | 怎么验证 |
|---|---|---|
| `adata.X` 已经是 log-normalized | ③+④ | 看到 `0.6931` = log1p(1) |
| `adata.layers['counts']` 原始 UMI | (③ 之前的备份) | 整数 vs X 是小数 |
| `adata.obs['n_counts'], n_genes` | ② QC | obs columns |
| `adata.obsm['X_pca']` | ⑥ PCA | 50 维 |
| `adata.obsp['distances'], connectivities` | ⑦ 表达 KNN | obsp |
| `adata.obs['Annotation']` | 已 annotate | 18 个 cell type |
| `adata.obs['spatial_leiden_e30_s8']` | 已 Leiden 聚类 | 30 NN, resolution 8 |
| `adata.obsp['spatial_connectivities']` | ⑨ 物理 KNN | obsp |
| `adata.uns['Annotation_colors']` | 配色已定 | 18 种颜色 |

所以拿到 ARTISTA 时，标准流程的 ①~⑨ 都不用自己跑。要做的是：
- **检查字段是否符合预期** (inspect)
- **可视化验证** (QC distributions, spatial scatter, UMAP)
- **加上 stGraphFM 特有的步骤** (构造 time label + 跨切片基因取交集)

---

## Raw `.h5ad` 完整流程代码骨架

以后看到没处理过的原始数据，套这个模板：

```python
import scanpy as sc
import squidpy as sq

# ① Inspect
adata = sc.read_h5ad('raw.h5ad')
print(adata)
print(adata.X[:5, :5].toarray())          # 是 int 还是 float？

# ② QC
adata.var['mt'] = adata.var_names.str.startswith('MT-')   # 哺乳动物
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.pct_counts_mt < 20]               # 过滤死细胞

# 备份原始 counts (NB / ZINB decoder 要用)
adata.layers['counts'] = adata.X.copy()

# ③ Normalize
sc.pp.normalize_total(adata, target_sum=1e4)

# ④ Log1p
sc.pp.log1p(adata)

# ⑤ HVG
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
adata.raw = adata                                          # 备份全部基因
adata = adata[:, adata.var.highly_variable]               # 只留 HVG

# ⑥ PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

# ⑦ Neighbors (表达空间)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# ⑧ UMAP + Leiden
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0)

# ⑨ Spatial neighbors (ST 特有)
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

# ⑩ Spatial autocorrelation (可选)
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100)
```

跑完后 `adata` 就是 ARTISTA 那种"开箱可用"的状态。

---

## 常见陷阱

| 陷阱 | 后果 | 处理 |
|---|---|---|
| 跳过 QC 直接 PCA | PC1 ≈ 测序深度 / 死细胞比例 | 先 ②③ 再 ⑥ |
| 没备份 `layers['counts']` | NB/ZINB decoder 没法训 | normalize 前先 `.copy()` 一份 |
| 直接在 `adata.X`（log 后）上算 Moran's I | 数值偏，但通常还能用 | 真要严谨：传 raw counts |
| 物种没有 `MT-` 前缀 | `pct_counts_mt` 全 0 | 蝾螈/水稻等：改 `qc_vars` 或跳过 mt 过滤 |
| 跨切片直接 concat | 各切片 HVG 不同 → `var_names` 不一致 | `anndata.concat(..., join='inner')` |
| 用 `obsp['connectivities']` 当 spatial graph | 那是表达 KNN，不是空间 KNN | 用 `spatial_connectivities` |

---

## 一句话总结

> 完整流程是 **inspect → QC → normalize → log1p → HVG → PCA → neighbors → UMAP**，ST 还要再加 **spatial neighbors + Moran's I**。**ARTISTA 已经跑完前 9 步**，日常使用只需要 inspect + 可视化。拿到没处理过的 raw `.h5ad` 时，从 ② QC 开始一步步来。
