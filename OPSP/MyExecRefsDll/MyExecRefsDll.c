// MyExecRefsDll.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
//#pragma warning(disable:4996)
int main()
{
	const char * location = "D:\\FangCloudSync\\Dropbox\\Github\\OPSP\\OPSP\\Debug\\";
	char * argv[1];
	argv[0] = (char*)location;
	wpsmain(1, argv);
    return 0;
}

