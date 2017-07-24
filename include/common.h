#pragma once
#ifndef _COMMON_H_
#define _COMMON_H_

#include<stdio.h>
#include<stdlib.h>
#define PAUSE system("pause");

//����������е���Сֵ
template<typename T>
T min(T a, T b, T c) {
	T temp = a < b ? a : b;
	return temp < c ? temp : c;
}

template<typename T>
T max(T  a, T  b) {
	return a > b ? a : b;
}

template<typename T>
T min(T a, T b) {
	return a < b ? a : b;
}


//ʹ��ģ�壬ָ�룬��̬�����ά���飬������ʱӦ��Ҫ�Ѿ���ĵ�ַ���������� &B
//template<typename T>
//void  allocMatrix(T*** pointer, int m, int n)
//{
//	(*pointer) = (T **)malloc((m + 1) * sizeof(T *));
//	if ((*pointer) == nullptr) {
//		printf("�ڴ�ռ䲻�㣬�޷������ά���飡\n");
//		exit(1);
//	}
//	for (int i = 0; i <= m; ++i) {
//		(*pointer)[i] = (T*)malloc((n + 1) * sizeof(T));
//		if ((*pointer)[i] == nullptr) {
//			printf("�ڴ�ռ䲻�㣬�޷������ά���飡\n");
//			PAUSE;
//			exit(1);
//		}
//	}
//}

//ʹ��ģ�壬ָ�붯̬�����ά���飬������ʱʹ�������û��ƣ�ֱ�ӽ����������������ɡ�
//���ַ�����������һ��ʹ�ô�ָ��ķ�����������
template<typename T>
void  allocMatrix(T**  &pointer, int m, int n)
{
	pointer = (T **)malloc((m + 1) * sizeof(T *));
	if (pointer == nullptr) {
		printf("�ڴ�ռ䲻�㣬�޷������ά���飡\n");
		exit(1);
	}
	for (int i = 0; i <= m; ++i) {
		pointer[i] = (T*)malloc((n + 1) * sizeof(T));
		if (pointer[i] == nullptr) {
			printf("�ڴ�ռ䲻�㣬�޷������ά���飡\n");
			PAUSE;
			exit(1);
		}
	}
}


//ʹ��ģ�壬ָ�붯̬�ͷŶ�ά���飬������ʱʹ�������û��ƣ�ֱ�ӽ����������������ɡ�
template<typename T>
void deleteMatrix(T ** &p, int m) {
	for (int i = m; i >= 0; --i) {
		free(p[i]);
		p[i] = nullptr;
	}
	free(p);
	p = nullptr;
}

//ʹ��ģ�壬ָ�붯̬������ά���飬������ʱʹ�������û��ƣ�ֱ�ӽ����������������ɡ�
template<typename T>
void  allocMatrix(T***  &pointer, int m, int n, int d)
{
	pointer = (T ***)malloc((m + 1) * sizeof(T **));
	if (pointer == nullptr) {
		printf("�ڴ�ռ䲻�㣬�޷������ά���飡\n");
		PAUSE;
		exit(1);
	}
	for (int i = 0; i <= m; ++i) {
		pointer[i] = (T**)malloc((n + 1) * sizeof(T*));
		if (pointer[i] == nullptr) {
			printf("�ڴ�ռ䲻�㣬�޷������ά���飡\n");
			PAUSE;
			exit(1);
		}
	}

	for (int i = 0; i <= m; ++i) {
		for (int j = 0; j <= n; ++j) {
			pointer[i][j] = (T*) malloc((d + 1) * sizeof(T));
			if (pointer[i][j] == nullptr) {
				printf("�ڴ�ռ䲻�㣬�޷������ά���飡\n");
				PAUSE;
				exit(1);
			}
		}
	}
}

//ʹ��ģ�壬ָ�붯̬�ͷŶ�ά���飬������ʱʹ�������û��ƣ�ֱ�ӽ����������������ɡ�
template<typename T>
void deleteMatrix(T *** &p, int m, int n) {
	for (int i = 0; i <= m; ++i) {
		for (int j = 0; j <= n; ++j) {
			free(p[i][j]);
			p[i][j] = nullptr;
		}
	}
	for (int i = m; i >= 0; --i) {
		free(p[i]);
		p[i] = nullptr;
	}
	free(p);
	p = nullptr;
}






#endif // !_COMMON_H_


