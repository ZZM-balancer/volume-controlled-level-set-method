// 定义参数
h_max = 0.04; // 最大网格尺寸
r = 0.1;     // 圆的半径

// 关键设置：控制保存哪些元素
Mesh.SaveAll = 0;          // 保存所有实体
Mesh.SaveParametric = 0;   // 不保存参数化信息
Mesh.SaveTopology = 0;     // 不保存拓扑信息
Mesh.MshFileVersion = 2.2;

// 定义立方体顶点 (1×2×1 立方体，x:0-1, y:0-2, z:0-1)
Point(1) = {0, 0, 0, h_max};
Point(2) = {1, 0, 0, h_max};
Point(3) = {1, 2, 0, h_max};
Point(4) = {0, 2, 0, h_max};
Point(5) = {0, 0, 1, h_max};
Point(6) = {1, 0, 1, h_max};
Point(7) = {1, 2, 1, h_max};
Point(8) = {0, 2, 1, h_max};

// 定义y=2平面上的圆
Point(9) = {0.5, 2, 0.5, h_max};    // 圆心
Point(10) = {0.5+r, 2, 0.5, h_max}; // 圆右侧
Point(11) = {0.5, 2, 0.5+r, h_max}; // 圆上侧
Point(12) = {0.5-r, 2, 0.5, h_max}; // 圆左侧
Point(13) = {0.5, 2, 0.5-r, h_max}; // 圆下侧

// 创建立方体的边
Line(1) = {1, 2};  // 底面边
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};  // 顶面边
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};  // 垂直边
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// 创建y=2平面上的圆边界
Circle(13) = {10, 9, 11}; // 右上弧
Circle(14) = {11, 9, 12}; // 左上弧
Circle(15) = {12, 9, 13}; // 左下弧
Circle(16) = {13, 9, 10}; // 右下弧

// 创建立方体的侧面(非y=2面)
Line Loop(1) = {1, 2, 3, 4};       // z=0底面
Plane Surface(1) = {1};

Line Loop(2) = {5, 6, 7, 8};       // z=1上面
Plane Surface(2) = {2};

Line Loop(3) = {9, -8, -12, 4};    // x=0后面
Plane Surface(3) = {3};

Line Loop(4) = {10, 6, -11, -2};   // x=1前面
Plane Surface(4) = {4};

Line Loop(5) = {1, 10, -5, -9};   // y=0左面
Plane Surface(5) = {5};

// 创建y=2平面上的两个表面
// 1. 圆内部表面
Line Loop(6) = {13, 14, 15, 16};   // 逆时针定义内部区域
Plane Surface(6) = {6};            // 圆内部

// 2. 圆外部表面 (y=2平面上挖去圆的区域)
Line Loop(7) = {-3, 11, 7, -12};       // 外部大矩形 (顺时针)
Line Loop(8) = {-13, -16, -15, -14}; // 内部圆边界 (顺时针，与外环方向相反)
Plane Surface(7) = {7, 8};         // 挖去圆的区域

// 创建体积，包含所有7个表面
Surface Loop(1) = {5, 1, 2, 3, 4, 7, 6}; // 6个侧面 + 2个顶面部分
Volume(1) = {1};

// 仅为面和体创建物理组（关键步骤）
Physical Surface(1) = {5};        // y=0
Physical Surface(2) = {6};        // y=2上的圆
Physical Surface(3) = {1, 2, 3, 4, 7}; // 其余所有表面
Physical Volume(4) = {1};         // 体积

// 生成3D网格
Mesh 3;

// 设置MSH文件版本

Save "Th008_010.msh";
Save "Th008_010.vtk";