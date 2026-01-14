面向OIer的NTT工具箱：NTT_Tools。

目前包含基础NTT、多项式乘法、多项式逆元、多项式开方、多项式求Ln、多项式求Exp。

可通过宏进行条件编译，**请注意不是一个功能对应一个宏，例如求Exp需要开启NTT_MUL、NTT_INV、NTT_LN、NTT_EXP四个宏**。

封装到namespace中，不会出现奇奇怪怪的冲突问题。