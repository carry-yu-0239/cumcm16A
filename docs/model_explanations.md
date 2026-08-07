# 模型说明索引

在本文件中汇总各子问题已经确认的模型说明，并链接至相应的模型卡、代码、输出和论文位置。

## 子问题一：静水风载下的准静态平衡

工作模型、参数口径和推导见 `model_cards/q1.md`。它以浮标吃水 \(\delta\) 为参数，
通过杆件静力递推和分段悬链线得到 \(H(\delta;v)\) 与 \(R(\delta;v)\)，再对
\(H=18\) m 求根。可复现实现及生成结果见 `src/q1/q1_static_equilibrium.py` 和
`outputs/q1/`；其中钢管浮力采用条件性口径，详见 `docs/open_questions.md` 的 OQ-1。
