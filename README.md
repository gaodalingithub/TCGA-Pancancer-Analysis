# TCGA-Pancancer-Analysis

### 项目概述
本流程是一个完整的TCGA（The Cancer Genome Atlas）泛癌单基因分析系统，用于系统性地探索目标基因在多种癌症类型中的表达特征、临床预后价值、分子机制及免疫治疗相关性。目标基因为 VCAN 。

| 序号 | 目录            | 功能描述                     |
| ---- | --------------- | ---------------------------- |
| 01   | 01.简介         | 项目介绍                     |
| 02   | 02.mRNAdownload | mRNA数据下载                 |
| 03   | 03.mRNAmerge    | mRNA数据合并                 |
| 04   | 04.idTrans      | ID转换（Ensemble ID转Symbol） |
| 05   | 05.clinical     | 临床数据下载                 |
| 06   | 06.getClinical  | 临床数据提取                 |
| 07   | 07.panDiff      | 泛癌差异分析                 |
| 08   | 08.diff         | 单基因差异分析               |
| 09   | 09.pairDiff     | 配对样本差异分析             |
| 10   | 10.survival     | 总生存期(OS)分析             |
| 11   | 11.PFS          | 无进展生存期(PFS)分析        |
| 12   | 12.ROC          | ROC曲线分析                  |
| 13   | 13.cliCor       | 临床相关性分析               |
| 14   | 14.cliHeatmap   | 临床特征热图                 |
| 15   | 15.Nomo         | 列线图(Nomogram)构建         |
| 16   | 16.indep        | 独立预后分析                 |
| 17   | 17.cor          | 基因相关性分析               |
| 18   | 18.circos       | Circos圈图可视化             |
| 19   | 19.diff         | 免疫细胞差异分析             |
| 20   | 20.GO           | GO功能富集分析               |
| 21   | 21.KEGG         | KEGG通路富集分析             |
| 22   | 22.GSEA         | GSEA基因集富集分析           |
| 23   | 23.estimate     | 肿瘤微环境评估               |
| 24   | 24.TMEvioplot   | 肿瘤微环境小提琴图           |
| 25   | 25.CIBERSORT    | 免疫细胞浸润分析             |
| 26   | 26.immuneCor    | 免疫相关性分析               |
| 27   | 27.Lollipop     | 棒棒糖图（体细胞突变）       |
| 28   | 28.checkpoint   | 免疫检查点相关性分析         |
| 29   | 29.TMBcor       | 肿瘤突变负荷相关性分析       |
| 30   | 30.pRRophetic   | 药物敏感性预测               |
| 31   | 31.TCIA         | TCIA肿瘤免疫评估             |
| 32   | 32.HPA          | HPA蛋白表达分析              |