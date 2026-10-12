// 定义参数
h_max = 0.06; // 最大网格尺寸

// 网格保存设置
Mesh.SaveAll = 0;
Mesh.SaveParametric = 0;
Mesh.SaveTopology = 0;
Mesh.MshFileVersion = 2.2;

// 定义长方体顶点
// 区域 D = (0,1) x (0,2) x (0,1)
Point(1) = {0, 0, 0, h_max};
Point(2) = {1, 0, 0, h_max};
Point(3) = {1, 2, 0, h_max};
Point(4) = {0, 2, 0, h_max};

Point(5) = {0, 0, 1, h_max};
Point(6) = {1, 0, 1, h_max};
Point(7) = {1, 2, 1, h_max};
Point(8) = {0, 2, 1, h_max};

// 定义长方体的边
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 1};

Line(5)  = {5, 6};
Line(6)  = {6, 7};
Line(7)  = {7, 8};
Line(8)  = {8, 5};

Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// z = 0 面
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// z = 1 面
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

// x = 0 面
Line Loop(3) = {9, -8, -12, 4};
Plane Surface(3) = {3};

// x = 1 面
Line Loop(4) = {10, 6, -11, -2};
Plane Surface(4) = {4};

// y = 0 面
Line Loop(5) = {1, 10, -5, -9};
Plane Surface(5) = {5};

// y = 2 面
Line Loop(6) = {3, 12, -7, -11};
Plane Surface(6) = {6};

// 创建体积
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// 物理边界编号
Physical Surface(1) = {5};          // y = 0
Physical Surface(2) = {6};          // y = 2
Physical Surface(3) = {1, 2, 3, 4}; // 其余表面

// 物理体编号
Physical Volume(4) = {1};

// 生成3D网格
Mesh.CharacteristicLengthMax = h_max;
Mesh.CharacteristicLengthMin = h_max;
Mesh 3;

// 保存文件
Save "Th006_box.msh";
Save "Th006_box.vtk";