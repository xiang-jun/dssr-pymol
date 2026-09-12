# DSSR-PyMOL GUI v1.3.3 中文快速说明

这是一个单文件 PyMOL 插件，使用 X3DNA-DSSR 检测 RNA/DNA 结构特征，并提供
图形界面、纯 Python RNA 二维编辑器、属性查询、JSON 查看，以及 2D/3D 双向选择。
主 GUI 和 RNA 2D studio 均使用统一的白色科学界面；开启果冻模式时背景仍为白色。

## 安装

1. 解压 ZIP。
2. 单独安装 `x3dna-dssr`。DSSR 没有包含在本压缩包中，请从 Columbia
   Technology Ventures 获取合法授权：
   <https://inventions.techventures.columbia.edu/technologies/CU20391>
3. 打开 PyMOL，进入 **Plugin → Plugin Manager → Install New Plugin**。
4. 选择本目录中的 `dssr_select.py`。
5. 完全关闭并重新启动 PyMOL。

如果不想永久安装，也可以在 PyMOL 命令行临时运行：

```text
run /完整路径/dssr_select.py
```

## 测试

```text
fetch 1ehz, async=0
dssr_gui
```

二维窗口既可以从主 GUI 打开，也可以直接运行：

```text
dssr_2d selection=1ehz
```

如果 PyMOL 找不到 DSSR，请在 GUI 的 **exe** 输入框填写完整路径，或者运行：

```text
dssr_2d selection=1ehz, exe=/完整路径/x3dna-dssr
```

Windows 路径建议使用正斜杠，例如：

```text
C:/Tools/DSSR/x3dna-dssr.exe
```

## 说明

- 2D 功能已经是纯 Python/Qt 实现，不需要安装 Jmol 或 Java。
- 发送给他人时，应发送完整 ZIP，不要单独发送 DSSR 可执行文件。
- 若用于论文或报告，请引用 `README.md` 中列出的插件与 DSSR 文献。
