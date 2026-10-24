"""Write evidence-bounded Chinese and English reports from computed metrics."""


def write_report(output,m):
    s=m['system']['case_stats']['interaction'];b=m['pressure_to_load']['delta'];p=m['paired_normal_budget'];a=m['attribution'];g=m['body_geometry']
    f=p['force_stats'];late=p['late_force_means']
    report=f'''# Orca v14：从压力差到载荷、再到机制的完整核验

**本轮结论：已有数据支持“体鳍组合改变上下表面压差，并产生正的系统载荷对照差”；还不支持“RV+ 单独造成全部增益”或“已定量证明涡量反转的来源”。流程已执行，不把物理预算未闭合掩盖为代码通过。**

数据范围是已存储的单个周期：B（Sweep35）9–10，I（Sweep35_OnlyFin）3–4，N（NoFin）9–10。三个算例各有 481 个载荷样本，压力/速度场预算采用 17 个相位。这不是多周期统计，也不是网格或时间收敛证明。

## 1. 先确定可比较的物理量

在当前 Neumann 压力设置下，每个算例允许独立的空间常数偏移 p_k(x,t)→p_k(x,t)+c_k(t)。所以采用 J_k=Cp_lower,k−Cp_upper,k，再计算 ΔJ=J_B−J_I。下表面在共同材料坐标的 x–z 图上做三角插值，不外推；最大边长 0.02L，固定配对随材料点运动。它是明确的几何配对定义，不冒充唯一的法向厚度配对。

共同区域是开口区域，直接积分其绝对压力仍然依赖零点。这里每个算例以共同下表面的面积加权平均压力 p_ref,k 作参考，并同时保留其区域面积向量。代码对场和参考一起施加独立、随相位变化的常数偏移，核验区域积分不变。闭合整体压力积分本身不受常数偏移影响；导出网格/接口的离散闭合误差仍需核验。

本报告使用向内的三角面积向量 A，以及 Cp=2p/U²、U=0.29、S_ref=0.548382。求解器力文件已经带两倍物理力的写出约定，故原始力列除以 U²S_ref，不能再多乘或少乘 2。这里的 p 是求解器运动学压力。

## 2. 效果：完整系统的对照差确实为正

定义 D_sys=C_L(B)−C_L(N)−C_L(I)，全部取闭合物体的载荷。

| 481 点周期平均 | D_sys |
|---|---:|
| 压力 | {s['pressure']['mean']:+.9f} |
| 黏性 | {s['viscous']['mean']:+.9f} |
| 文件中的总力 | {s['total']['mean']:+.9f} |

压力正增量大于黏性负增量，因此该周期总的对照差为正。文件分别四舍五入压力、黏性和总力，三者不要求逐点机器精度相加；最大差为 {m['system']['writer_sum_max_error']:.3g}。

17 相位的压力平均值为 {m['system']['sampling']['pressure']['sparse_mean']:.9f}，比 481 点值高 {100*m['system']['sampling']['pressure']['sparse_mean_error']/s['pressure']['mean']:.2f}%。后文 17 相位局部预算不能冒充高时间分辨率的全周期结论。

B/N 的共享顶点运动也重新读取了原始标记面：以 5×10⁻⁶L 的初始匹配容差匹配到 {g['matched_vertices']}/{g['nofin_vertices']} 个顶点，匹配部分全周期最大位置差 {g['max_matching_motion_error']:.3g}L。未匹配顶点可能来自不同剖分或不同几何；这一检查没有证明全体表面运动完全相同。因此 D_sys 是这三个已存储算例的组合对照差，对“仅改变体鳍耦合”的因果解释仍以完整运动/几何控制核验为条件。

![整体对照与配对压力载荷](09_system_and_pair.png)

## 3. 桥梁：上下压差确实能解释同一区域的载荷变化

定义 F_pair=Σ ΔJ·(−A_y,upper)/S_ref。两个鳍共 {m['pressure_to_load']['query_faces']} 个上表面查询三角面。该量是配对压差的薄面载荷代理，有限厚度、未配对部分及几何插值不能悄悄略去。

| 共同区域、同一参考下的周期平均变化 | 数值 |
|---|---:|
| 配对压差积分 F_pair | {b['stats']['pair']['mean']:+.9f} |
| 三角标记面的参考压力积分 | {b['stats']['marker_reference']['mean']:+.9f} |
| 原始笛卡尔接口的参考压力积分 | {b['stats']['native_reference']['mean']:+.9f} |
| 几何、未配对与插值余项 | {b['stats']['geometry_and_unpaired']['mean']:+.9f} |
| 笛卡尔接口减标记面积分 | {b['stats']['discretization']['mean']:+.9f} |

逐相位核验：接口参考力=配对压差积分+几何/未配对/插值余项+离散化差。配对积分相对接口参考力的 RMS 误差为其波动标准差的 {100*b['pair_vs_native']['rms_over_measured_fluctuation']:.2f}%；标记面积分相对接口仅为 {100*b['marker_vs_native']['rms_over_measured_fluctuation']:.2f}%。这是有用的压力—载荷桥梁，但它仍然是指定参考下的开口区域载荷。

重新读取的 B/I 全部 34 个原始压力接口，映射后与已存储缓存逐项一致。该一致性验证读取和映射复现，不是独立验证求解器内部边界处理。完整接口重建相对实际求解器压力力的波动归一化 RMS 误差分别为 {100*m['pressure_to_load']['baseline']['whole_reconstruction']['rms_over_measured_fluctuation']:.3f}% 和 {100*m['pressure_to_load']['isolated']['whole_reconstruction']['rms_over_measured_fluctuation']:.3f}%；内部六向边界掩码、物体归属和双精度几何量未全部导出。

**不能跨越的推论：** F_pair 的平均值接近 D_sys 的平均值，不代表它解释了全部系统力。两者在同样 17 相位上相关系数为 {m['pair_vs_system']['correlation']:.4f}，但 RMS 差达到系统波动标准差的 {100*m['pair_vs_system']['rms_over_measured_fluctuation']:.1f}%。区域范围不同、身体响应及未配对部分仍在其中；这个差不是单一数值误差，也不能直接命名为 RV 力。

![同一区域的压力到载荷桥梁](10_pressure_to_load.png)

## 4. 压差从何而来：保留外层、解析层和壁面间隙

负 z 一侧的鳍有 {p['faces']} 个共同有效配对点，覆盖原共同上表面面积的 {100*p['common_upper_area_coverage']:.2f}%。对上下表面分别沿外法向采样，再做配对相减。逐面、逐相位成立：

J_surface = J_outer(h₁) + [J(h₀)−J(h₁)] + [J_surface−J(h₀)]，其中 h₀=0.012L，h₁=0.040L。

| 单鳍配对载荷变化 | 周期均值 | 后四分之一周期均值 |
|---|---:|---:|
| 表面总压差 | {f['surface']['mean']:+.9f} | {late['surface']:+.9f} |
| 外层配对压力 | {f['outer']['mean']:+.9f} | {late['outer']:+.9f} |
| h₀–h₁ 解析层压差 | {f['layer']['mean']:+.9f} | {late['layer']:+.9f} |
| 壁面–h₀ 间隙 | {f['wall']['mean']:+.9f} | {late['wall']:+.9f} |

这是压力的空间记账恒等式，不是三个独立因果贡献。尤其不能把外层压力称为“没有涡的背景”；外层压力同样由整个流场决定。

不可压常黏度动量方程给出 −∇p=∂u/∂t+∇(|u|²/2)−u×ω−ν∇²u。沿两侧法向分别积分、再配对相减，得到解析层的 Lamb、动能梯度、非定常和黏性项。移动采样射线的时间导数已经减去 (∇u)·v_ray；未进行拟合、幅值缩放或以余项校准各项。

| 解析层预算 | 周期均值 | 后四分之一周期均值 |
|---|---:|---:|
| Lamb | {f['lamb']['mean']:+.9f} | {late['lamb']:+.9f} |
| 动能梯度 | {f['kinetic']['mean']:+.9f} | {late['kinetic']:+.9f} |
| 非定常 | {f['unsteady']['mean']:+.9f} | {late['unsteady']:+.9f} |
| 黏性 | {f['viscous']['mean']:+.9f} | {late['viscous']:+.9f} |
| 四项合计 | {f['prediction']['mean']:+.9f} | {late['prediction']:+.9f} |
| 实测解析层 | {f['layer']['mean']:+.9f} | {late['layer']:+.9f} |
| 实测减重建 | {f['closure']['mean']:+.9f} | {late['closure']:+.9f} |

正 Lamb 项被动能和非定常负项显著抵消，后四分之一周期尤其明显。因此不能用 Lamb 项大小直接代表净压差或某个涡的升力。周期平均误差为实测层均值的 {100*f['closure']['mean']/f['layer']['mean']:.2f}%，但瞬态 RMS 误差仍为波动标准差的 {100*p['layer_reconstruction']['rms_over_measured_fluctuation']:.1f}%。换用中心时间差分后该误差为 {100*m['central_derivative_sensitivity']['layer_reconstruction']['rms_over_measured_fluctuation']:.1f}%，提示当前 17 相位导数敏感；两种方法的差异不是时间收敛证据。ν 来自现有输入文件，运行副本身份尚未独立归档。

![配对压力预算与抵消项](11_paired_pressure_budget.png)

## 5. 涡结构证据能走到哪一步

前轮六条完整的配对轨迹支持 B 中 ω_x 由负转正、I 中仍为负；这些相关轨迹不构成六个独立重复样本，也不证明整个涡管翻转。RV+ 的环量与核–环压力差支持局部旋转及相对低压的共同出现。但前轮 Neumann 一致的空间对照显示，上游等面积区域的后期压差增量 0.59573，高于 RV+ 邻域的 0.45163。因此 RV+ 邻域不是唯一增强位置，当前证据更符合分布式压力重分配，而非独占的局部吸力来源。该判断是对现有空间对照的解释，未建立唯一因果机制。

本轮不是只说“需要力分解”：已经重新审计 {a['saved_potential_results_checked']} 个保存的辅助势结果。势求解最大相对残差 {a['max_potential_solver_residual']:.3g}，目标记账误差 {a['max_target_identity_error']:.3g}，但宽域、去两层近壁单元方案仍有未分配均值 {a['provisional_delta']['unassigned']:+.9f}，相当于其原区域目标的 {100*a['unassigned_over_target']:.1f}%。其中约 93.35% 是定义上的近壁余项。原目标为单鳍区域、旧参考约定，不能把这些比例移植到本轮闭合系统 D_sys。数值方程解得准不等于物理归因闭合。

涡量来源的 {a['transport']['windows']} 个已算窗口中，仅 {a['transport']['all_terms_valid_windows']} 个全部项采样有效；即使只看有效窗口，绝对余项/|Δω_x| 的中位数仍为 {100*a['transport']['valid_window_median_abs_residual_over_change']:.1f}%。所以暂不能认定拉伸、倾斜或黏性中的哪一项已经定量解释反号。

## 6. 论文应使用的逻辑，以及下一项真正要解决的问题

可写的主线是：**三算例的闭合载荷对照 → Neumann 不变的上下压差 → 同一区域的压力积分载荷 → 外层/解析层/壁面间隙 → 解析层内部相互抵消的动量项 → 涡结构的空间与时间证据 → 明确未闭合的来源与单涡归因。** 这条链已经有结果、公式、三张新图、输入哈希和一键运行代码；最后两步不能用一张涡图补齐。

文献的作用已经落实到分析设计：[组内及方法文献矩阵](../2026-10-23-mechanism-logic/LITERATURE_MATRIX.md) 中 L01 对应三模型对照，L05 对应“效果→相位→结构→压力→力”的英文组织，L08 对应同时核验均值和波形，L07 对应保留应变和边界项的力分解。本轮没有新增全文阅读，不把先前选定段落阅读夸大为全部文献已审查。

优先后续工作是将压力与速度导数预算改为求解器原生时间差分/更密相位输出，并导出近壁边界掩码与运动信息；同一流程直接检验余项是否下降。现有三种辅助势区域已经算过，继续画同类 Q 图无法消除近壁余项。正式单涡归因之前，还要将区域目标转换为一致参考或闭合控制面，并做域、近壁剔除和时空分辨率收敛。现有结果没有证明新的高分辨率计算已经执行。

## 运行与交付

运行 `python code/run__researchPipeline.py`，依次执行 34 相位原始接口重读、481 点载荷对照、压差载荷桥梁、配对动量预算、几何检查、102 个既有势解与涡量窗口审计、PNG/PDF/SVG 出图和本报告生成。配置路径位于两个被调用运行脚本的 main()。依赖 numpy/scipy/matplotlib/pyvista，以及本地 Picar_Post 的 TecIO 读取器和 DLL；原始数据留在本地，不上传 GitHub。`chain_signals.npz` 是本轮导出的少量报告曲线，不是原始 CFD 数据。

结果：[完整指标](integrated_chain.json)、[原始接口核验](raw_interface_audit.json)、[流程执行记录](pipeline_run.json)、[英文结果段落](MANUSCRIPT_DRAFT.md)、[已执行清单](TODO.md)。v12 单侧绝对压力的“19.3% 更强吸力”继续撤回，不能引用。
'''
    (output/'README.md').write_text(report,encoding='utf-8')
    english=f'''# Results draft — pressure differences, loads and bounded attribution

For the stored cycle, the closed-body contrast B−N−I yielded mean pressure, viscous and total lift coefficients of {s['pressure']['mean']:.7f}, {s['viscous']['mean']:.7f} and {s['total']['mean']:.7f}, respectively. The pressure contrast remained positive after the viscous contribution was included. These values use 481 force samples per case; using only the 17 field-output phases increased the mean pressure contrast by 9.95%. This comparison characterizes the three stored simulations. Shared-vertex motion was checked, but complete surface-kinematic equivalence and multi-cycle convergence have not been established.

To remove the independent pressure constants permitted by the Neumann formulation, we compared paired lower-to-upper pressure differences on the common material surface. Their integrated lift proxy had a mean B−I increment of {b['stats']['pair']['mean']:.7f}. With an explicit lower-surface pressure reference, the corresponding Cartesian-interface regional load increment was {b['stats']['native_reference']['mean']:.7f}. The paired proxy differed from that regional signal by 8.96% of its fluctuation standard deviation; including the finite-geometry, unpaired-surface and interpolation remainder reduced the marker-to-interface discrepancy to 1.03%. In contrast, comparison against the whole-system signal produced a 75.1% normalized RMS difference despite a correlation of 0.9926. Agreement of the cycle averages therefore does not establish that the common-fin pressure increment accounts for the whole-system response.

For one fin, the jointly valid pressure pairs covered 96.24% of the common upper-surface area. The paired surface-pressure increment was separated exactly into an outer-pressure term, a resolved normal-layer difference and a wall-gap term. Their mean lift-proxy contributions were {f['outer']['mean']:.7f}, {f['layer']['mean']:.7f} and {f['wall']['mean']:.7f}. This is spatial pressure bookkeeping, not a causal partition. Within the resolved layer, the positive Lamb contribution was opposed by kinetic-gradient and unsteady terms. The reconstructed mean was within 3.18% of the measured mean, while the phase-resolved error remained 39.2% of the measured fluctuation standard deviation and increased to 61.2% with a central time derivative. These results support a balance among competing terms, rather than identification of the Lamb term with a vortex-specific lift contribution.

Previously computed auxiliary-potential results passed the algebraic target checks but retained an unassigned regional mean equal to −55.6% of their original target. That target and pressure convention differ from the present closed-system contrast. Only 44 of 98 stored vorticity-transport windows had valid samples for all terms; their median absolute residual relative to the vorticity change was 47.7%. Together with the stronger paired-pressure increment in an equal-area upstream control than near RV+, these findings support distributed pressure redistribution associated with body–fin interaction, but do not establish a unique RV+ force fraction or a quantitatively closed origin of the streamwise-vorticity sign change.

Provenance and limitations: this is a draft based on one stored cycle, not a validated manuscript conclusion. The local viscosity input has not been independently matched to the run copy. See the Chinese report for equations, pressure-reference conventions, figure captions and numerical checks. Literature reading scope is recorded in the existing literature matrix; no additional full-text reading is claimed here.
'''
    (output/'MANUSCRIPT_DRAFT.md').write_text(english,encoding='utf-8')
    (output/'TODO.md').write_text('''# Executed workflow and remaining scientific tests

- [x] Inventory force histories, surface caches, raw interfaces, profiles and readers.
- [x] Use Neumann-invariant paired pressure and explicitly referenced regional loads.
- [x] Re-read all 34 B/I raw pressure interfaces; compare each mapped cache.
- [x] Integrate 481-sample closed-system pressure, viscous and total contrasts; compare 17-sample means.
- [x] Connect paired pressure to regional marker and Cartesian loads; retain geometric and discretization remainders.
- [x] Compute paired upper/lower normal-momentum budgets and derivative sensitivity.
- [x] Read NoFin marker motion; report matched coverage instead of assuming full equivalence.
- [x] Audit 102 saved auxiliary potentials and all 98 vorticity windows.
- [x] Generate figures, Chinese report, English results and machine-readable metrics through one runner.
- [ ] Establish full geometry/motion equivalence and multi-cycle/grid/time convergence.
- [ ] Close the native near-wall momentum and vorticity budgets with independently archived run parameters.
- [ ] Rebuild a reference-consistent vortex-force target before claiming an RV+ force fraction.

Unchecked items are unresolved scientific requirements, not claims that the available analysis code was never run. No production CFD rerun or new elliptic solve is claimed by this pipeline.
''',encoding='utf-8')
