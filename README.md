# Orca — 代码与研究报告

本仓库保存助手为 Orca 编写或修改的代码及对应报告，供远程查看。仿真数据、个人文档和论文全文保留本地。

## 从这里开始

- **[v13：Neumann 边界下的压差关系与修正](reports/2026-10-23-pressure-gauge/README.md)**
- [压差机制英文修订](reports/2026-10-23-pressure-gauge/MANUSCRIPT_NOTE.md)

- **[v12 前版研究报告（压力解释已修正）：文献、空间对照与机制论证](reports/2026-10-23-mechanism-logic/README.md)**
- **[8 篇原文的证据矩阵与推荐阅读顺序](reports/2026-10-23-mechanism-logic/LITERATURE_MATRIX.md)**
- **[按组内论证方式重写的英文草稿](reports/2026-10-23-mechanism-logic/MANUSCRIPT_DRAFT.md)**
- [v11：平均/瞬态闭合、配对轨迹与局部压力](reports/2026-10-23-evidence-audit/README.md)
- [23 份去重历史报告目录](reports/history/README.md)：包含旧版、草稿和片段；最新结论优先。
- [下一步：压力到完整鳍载荷的闭合](reports/2026-10-23-mechanism-logic/NEXT.md)

本轮结论：以各算例内部的上下表面压差 J 为主变量，再比较 ΔJ。共同配对区域周期平均 J 为体鳍组合 0.41757、孤立鳍 0.26030；时变压力偏移检验保持不变。组合算例的鳍标签区域在根部开放，不能直接将其压力积分与闭合孤立鳍比较而忽略参考约定。

![上下表面压差](reports/2026-10-23-pressure-gauge/08_pressure_jump.png)

## 代码与复算

- `code/run__pressureJump.py`：材料投影配对、跨鳍压差、压力常数不变性和开放面片敏感性，图 8。

- `code/run__evidenceAudit.py`：复核均值、相位闭合和轨迹，图 1–2。
- `code/run__localPressureBridge.py`：固定材料邻域压力与动量预算，图 3–4。
- `code/run__spatialControls.py`：等面积空间对照、非重叠敏感性、压力分布和 Lamb 中点展开，图 5–7。
- `code/vortex_atlas/run__vortexAtlas.py`：v10 涡图维护代码快照；本轮未重跑其流程。

先在 runner 的 main() 设置本机输入路径。新分析依赖 NumPy、Matplotlib 和 SciPy。原始输入未上传；JSON 中记录输入哈希，结果不代表无输入数据也能复算。

图内数据以栅格层嵌入 PDF/SVG，文字和坐标保留可编辑元素，另有 600 dpi PNG。GitHub 可直接阅读 Markdown 与 PNG；下载合集后可打开 INDEX.html 查看图件及历史报告。代码内部一致性检查通过不等于 CFD 网格收敛或实验验证。

项目规则见 [PROJECT.md](docs/agent-loop/PROJECT.md)。报告与英文草稿均非投稿定稿。
