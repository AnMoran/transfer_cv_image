# transfer_cv_image

1.本工程采用纯C代码实现了 opencv中的4点的透视变换，和基于5点的人脸归一化
    由于工程中通常需要把代码移植到前端，安卓上可以通过编译armv7或者armv8版本的opencv来使用opencv中的一些图像操作，但是
    IOS端由于不支持动态库的加载，而opencv.framework的体积又过于庞大，因此经常需要手动实现一些opencv的比较困难的图像操作
    比如3点放射变换归一化图像，由于三点的放射变换较为简单，本工程中只实现了4点的透视变换，经常用于身份证、车牌等文本行的归一化
    和5点的人脸归一化的纯c语言实现

2.基于4个坐标点的透视变换实现如下
    norm_4pt.h
    norm_4pt.cpp
    在main.cpp中采用opencv验证实现的是否正确
    验证发现求取的透视变换矩阵和opencv一样，但是归一化后的图像还稍微有差异，但是已不太影响。

3.基于5点的人脸归一化
    todo
