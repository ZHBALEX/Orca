# Orca — 代码与研究报告

本仓库保存助手为 Orca 编写或修改的代码及对应报告，供远程查看。仿真数据、个人文档和论文全文保留本地。

## 从这里开始

- **[v20：全周期演化与上侧流动来源回溯](reports/2026-10-24-dong-evolution/README.md)**

- **[v19：Dong 组路线的同相位流场—压差对照](reports/2026-10-23-dong-comparison/README.md)**

- **[v18：固定周期内的压差与增升机制](reports/2026-10-23-phase-mechanism/README.md)**

- **[v17：周期稳定性、采样误差与发表证据边界](reports/2026-10-23-cycle-evidence/README.md)**
- [v17 英文结果草稿](reports/2026-10-23-cycle-evidence/MANUSCRIPT_DRAFT.md)

- **[v16：把旋转相关压差增量追到剪切与法向速度变化](reports/2026-10-23-rotation-origin/README.md)**
- [v16 英文论证](reports/2026-10-23-rotation-origin/MANUSCRIPT_DRAFT.md)

v16 结果：在同一1071个上下配对面上，正Lamb增量由涡量变化项提供，速度变化项部分抵消；进一步定位到上侧解析层的法向剪切和法向速度弦向梯度。NoFin同位置压力向近表面减弱，而B−I增量增强。分析范围和动量余项见v16。


- **[v15：为什么平均鳍载荷增大——四算例对照与机制解释](reports/2026-10-23-lift-cause-controls/README.md)**
- [v15 英文论证](reports/2026-10-23-lift-cause-controls/MANUSCRIPT_DRAFT.md)

v15 结果：将运动冻结于平均姿态后，共同鳍面压差载荷的身体效应下降 76.55%；内侧 60% 共同展向贡献 84.73%。增益体现为平均压差改变，平均迎角增大与直接加速度均值均不能解释；已解析旋转相关压力源项在六种设置中保持正贡献。具体范围与限定见 v15。


- **[v14：压力差—载荷—机制的完整核验](reports/2026-10-23-integrated-chain/README.md)**
- [v14 英文结果](reports/2026-10-23-integrated-chain/MANUSCRIPT_DRAFT.md) · [实际执行记录](reports/2026-10-23-integrated-chain/pipeline_run.json) · [已执行清单与未闭合问题](reports/2026-10-23-integrated-chain/TODO.md)

v14 完整跑通原始接口重读、密集载荷、配对压力载荷、动量预算及既有力分解审计。系统平均压力对照差 +0.006120；区域压差能较好解释同范围载荷，但不能解释整个系统波形。正 Lamb 项受到动能和非定常项抵消，单涡因果归因仍未闭合。


- **[v13：Neumann 边界下的压差关系与修正](reports/2026-10-23-pressure-gauge/README.md)**
- [压差机制英文修订](reports/2026-10-23-pressure-gauge/MANUSCRIPT_NOTE.md)

- **[v12 前版研究报告（压力解释已修正）：文献、空间对照与机制论证](reports/2026-10-23-mechanism-logic/README.md)**
- **[8 篇原文的证据矩阵与推荐阅读顺序](reports/2026-10-23-mechanism-logic/LITERATURE_MATRIX.md)**
- **[按组内论证方式重写的英文草稿](reports/2026-10-23-mechanism-logic/MANUSCRIPT_DRAFT.md)**
- [v11：平均/瞬态闭合、配对轨迹与局部压力](reports/2026-10-23-evidence-audit/README.md)
- [23 份去重历史报告目录](reports/history/README.md)：包含旧版、草稿和片段；最新结论优先。
- [v12 历史下一步（最新进展见 v14）](reports/2026-10-23-mechanism-logic/NEXT.md)

v13 结论：以各算例内部的上下表面压差 J 为主变量，再比较 ΔJ。共同配对区域周期平均 J 为体鳍组合 0.41757、孤立鳍 0.26030；时变压力偏移检验保持不变。组合算例的鳍标签区域在根部开放，不能直接将其压力积分与闭合孤立鳍比较而忽略参考约定。

![上下表面压差](reports/2026-10-23-pressure-gauge/08_pressure_jump.png)

## 代码与复算

- `code/run__pairedRotationOrigin.py`：v16 中点/方向/梯度精确分解、共同有效面敏感性、外层至表面对照、图14–15；实际运行通过。


- `code/run__liftCauseControls.py`：v15 四算例运动/静止对照、压差与投影恒等式、来流检验、源项敏感性和图 12–13；已运行通过。


- `code/run__researchPipeline.py`：v14 一键流程，原始数据→压力/载荷→预算审计→图 9–11→中英文报告；本机完整执行通过。


- `code/run__pressureJump.py`：材料投影配对、跨鳍压差、压力常数不变性和开放面片敏感性，图 8。

- `code/run__evidenceAudit.py`：复核均值、相位闭合和轨迹，图 1–2。
- `code/run__localPressureBridge.py`：固定材料邻域压力与动量预算，图 3–4。
- `code/run__spatialControls.py`：等面积空间对照、非重叠敏感性、压力分布和 Lamb 中点展开，图 5–7。
- `code/vortex_atlas/run__vortexAtlas.py`：v10 涡图维护代码快照；本轮未重跑其流程。

先在 runner 的 main() 设置本机输入路径。新分析依赖 NumPy、Matplotlib、SciPy 和 PyVista；原始接口读取另需本机 Picar_Post/TecIO 与 DLL。原始输入未上传；JSON 中记录输入哈希，结果不代表无输入数据也能复算。

图内数据以栅格层嵌入 PDF/SVG，文字和坐标保留可编辑元素，另有 600 dpi PNG。GitHub 可直接阅读 Markdown 与 PNG；下载合集后可打开 INDEX.html 查看图件及历史报告。代码内部一致性检查通过不等于 CFD 网格收敛或实验验证。

项目规则见 [PROJECT.md](docs/agent-loop/PROJECT.md)。报告与英文草稿均非投稿定稿。
