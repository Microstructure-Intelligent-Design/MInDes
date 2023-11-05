/*
This file is a part of the microstructure intelligent design software project.

Created:     Qi Huang 2023.04

Modified:    Qi Huang 2023.04;

Copyright (c) 2019-2023 Science center for phase diagram, phase transition, material intelligent design and manufacture, Central South University, China

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free
	Software Foundation, either version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once
#include "../../solvers/base/sysTool.h"
using namespace std;
namespace pf {
	class BMP24reader
	{
	public:
		BMP24reader() {};
		~BMP24reader() {
			delete bmpbuf;
			bmpbuf = nullptr;
		};
		vector<int> getRGB(int x, int y) {
			vector<int> rgb;
			rgb.push_back(bmpbuf[y * linebyte + x * 3 + 2]);
			rgb.push_back(bmpbuf[y * linebyte + x * 3 + 1]);
			rgb.push_back(bmpbuf[y * linebyte + x * 3 + 0]);
			return rgb;
		}
		double getGrayPercentage(int x, int y) {
			vector<int> rgb;
			rgb.push_back(bmpbuf[y * linebyte + x * 3 + 2]);
			rgb.push_back(bmpbuf[y * linebyte + x * 3 + 1]);
			rgb.push_back(bmpbuf[y * linebyte + x * 3 + 0]);
			return (rgb[0] * 0.299 + rgb[1] * 0.587 + rgb[2] * 0.114) / 255;
		}
		void changeRGB(int x, int y, int R, int G, int B) {
			bmpbuf[y * linebyte + x * 3 + 2] = R;
			bmpbuf[y * linebyte + x * 3 + 1] = G;
			bmpbuf[y * linebyte + x * 3 + 0] = B;
		}
		void GrayTransfer() {
			for(int x = 0; x < bmp_width; x++)
				for (int y = 0; y < bmp_height; y++) {
					vector<int> rgb = getRGB(x, y);
					int gray = double_to_int(rgb[0] * 0.299 + rgb[1] * 0.587 + rgb[2] * 0.114);
					changeRGB(x, y, gray, gray, gray);
				}
		}
		//检查打开的文件路径是否正确，是否是bmp文件
		int safe(string fileName) 
		{
			const char* bmpname = fileName.c_str();
			BITMAPINFOHEADER head;  //BITMAPINFOHEADER为编译系统自带的结构体，故直接使用
			FILE* file = fopen(bmpname, "rb");//打开文件，并且以二进制方式读取
			if (!file)               //文件是否能打开
			{
				printf("file dont exist, please check !.\n");
				return 0;//0默认为打开的文件不正确
			}
			if (!(fgetc(file) == 'B' && fgetc(file) == 'M'))//打开的是否是bmp图片
			{
				printf("file isnt in bmp format, please check.\n");
				return 0;
			}//其余的都默认为正确的bmp文件
			fseek(file, sizeof(BITMAPFILEHEADER), 0);//跳过位图文件头
			fread(&head, sizeof(BITMAPINFOHEADER), 1, file);//将位图信息头读入head结构体中。
			if (head.biBitCount != 24)
			{
				printf("The program now only supports true color BMP files\n");
				return 0;
			}
			fclose(file);
			return 1;//进行到这一步，文件默认为24位bmp文件
		}
		//将数据读入
		int read(string fileName)
		{
			const char* bmpname = fileName.c_str();
			BITMAPINFOHEADER head;//head为位图信息头结构体
			FILE* file = fopen(bmpname, "rb");//以二进制方式打开
			fseek(file, sizeof(BITMAPFILEHEADER), 0);//跳过位图文件头
			fread(&head, sizeof(BITMAPINFOHEADER), 1, file);//将位图信息头读入head
			bmp_width = head.biWidth;//图片宽
			bmp_height = head.biHeight;//图片高
			bibitcount = head.biBitCount;//图片颜色表
			linebyte = (bmp_width * bibitcount + 31) / 32 * 4;//每行像素所占的字节数，
		//Windows规定一个扫描行所占的字节数必须是4的倍数（即以unsigned int为单位），不足的以0填充，0补在每一行的最后
		//因此位图数据所占用的字节总数并不是图像的高度乘以图像宽度再乘以每像素所占字节数就行了，正确的算法应该如前所示
			bmpbuf = (unsigned char*)malloc(linebyte * bmp_height);//开辟空间用以保存
			fread(bmpbuf, 1, linebyte * bmp_height, file);//读取位图的数据
			fclose(file);
			return 1;
		}
		//将修改后的图片写入另一个地方
		int save(string newfileName)
		{
			const char* newbmpname = newfileName.c_str();
			int colortablesize = 0;//颜色表大小，彩色图像为0
			BITMAPFILEHEADER filehead;//位图文件头
			BITMAPINFOHEADER datehead;//位图信息头
			FILE* file = fopen(newbmpname, "wb");//打开文件并以二进制方式写入
			if (!file) return 0;//打开不成功
			filehead.bfType = 0x4D42;//文件为bmp格式
			filehead.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + colortablesize + linebyte * bmp_height;//位图文件大小
			filehead.bfReserved1 = 0;//位图文件保留字，必须为0
			filehead.bfReserved2 = 0;//位图文件保留字，必须为0
			filehead.bfOffBits = 54 + colortablesize;//位图数据的起始位置
			fwrite(&filehead, sizeof(BITMAPFILEHEADER), 1, file);//写入位图文件头
			datehead.biBitCount = bibitcount;//像素所需位数
			datehead.biClrImportant = 0;//位图显示过程中重要的颜色，0表示全部都重要
			datehead.biClrUsed = 0;//颜色表中颜色数
			datehead.biCompression = 0;//压缩类型，必须是0
			datehead.biHeight = bmp_height;//位图高度
			datehead.biPlanes = 1;//目标设备级别
			datehead.biSize = 40;//该结构体字节数
			datehead.biSizeImage = linebyte * bmp_height;//位图大小
			datehead.biWidth = bmp_width;//位图宽度
			datehead.biXPelsPerMeter = 0;//位图水平分辨率
			datehead.biYPelsPerMeter = 0;//位图垂直分辨率
			fwrite(&datehead, sizeof(BITMAPINFOHEADER), 1, file);//写入位图信息头
			fwrite(bmpbuf, linebyte * bmp_height, 1, file);//写入位图数据
			fclose(file);//写入完毕
			return 1;
		}
		int bmp_width;               //宽
		int bmp_height;              //高
	private:
		/*下列变量经常使用，故写为全局变量*/
		int linebyte;               //计算每行所占的字节数
		int bibitcount;             //颜色表
		unsigned char* bmpbuf;      //位图图像的数据
	};
}