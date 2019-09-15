
disp('To mex random forest by C++');

mex src/mex_ClassificationRF_train.cpp  src/classRF.cpp src/classTree.cpp src/rfutils.cpp src/cokus.cpp src/buildtree.cpp -output mexClassRF_train  -lm -DMATLAB -O 
mex src/mex_ClassificationRF_predict.cpp src/classRF.cpp src/classTree.cpp src/rfutils.cpp src/cokus.cpp src/buildtree.cpp -output mexClassRF_predict -lm -DMATLAB -O

disp('Cool, mex done.');

