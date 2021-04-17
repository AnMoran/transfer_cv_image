#include"norm_5pt.h"
#include"common.h"
//using namespace std;

int NormImg_5_points(const float* pImg, int nImgWidth, int nImgHeight,int channels,
	std::vector<std::vector<float>> src,std::vector<std::vector<float>> dst,float *pfNormImg)
{
	/*********************************************/
	int i, j, k;
	int m = 2;
	int n = 5;

	std::vector<std::vector<float>> src_new;
	std::vector<std::vector<float>> src_matrix;	
	std::vector<std::vector<float>> dst_new;
	std::vector<std::vector<float>> dst_matrix;	

	center_and_normalize_points(src, src_new, src_matrix);
	center_and_normalize_points(dst, dst_new, dst_matrix);

	std::vector<std::vector<float> > A(n * 2, std::vector<float>(5, 0));
	for (i = 0; i < n; i++)
	{
		A[i][0] = src_new[i][0];
		A[i][1] = 1;
		A[i][2] = -src_new[i][1];
		A[i][4] = dst_new[i][0];
		A[i + n][2] = src_new[i][0];
		A[i + n][0] = src_new[i][1];
		A[i + n][3] = 1;
		A[i + n][4] = dst_new[i][1];
	}

	//svd
	std::vector<std::vector<float> > U;
	std::vector<float> S;
	std::vector<std::vector<float>> V;
	svd(A, 5, U, S, V);
	int nCnt = 0;
	while (fabs(V[4][4]) < 1e-3 && nCnt < 100)
	{
		nCnt ++;
		svd(A, 5, U, S, V);
	}


	std::vector<std::vector<float> >H(3, std::vector<float>(3,0));
	if (fabs(V[4][4]) < 1e-3)
		return -1;
	float a0 = -V[4][0] / V[4][4];
	float a1 = -V[4][1] / V[4][4];
	float b0 = -V[4][2] / V[4][4];
	float b1 = -V[4][3] / V[4][4];
	H[0][0] = a0;
	H[0][1] = -b0;
	H[0][2] = a1;
	H[1][0] = b0;
	H[1][1] = a0;
	H[1][2] = b1;
	H[2][2] = 1;

	std::vector<std::vector<float>> AA;
	if (Gauss(dst_matrix, AA) == false)
	{
		return -1;
	}

	
	std::vector<std::vector<float>> C; 
	std::vector<std::vector<float>> D;
	CalMatDot(H, src_matrix, C);
	CalMatDot(AA, C, D);

	//printf("D\n");
	print(D);

	warpAffineLinear(pImg, pfNormImg, nImgWidth, nImgHeight, NORM_WIDTH, NORM_HEIGHT, channels, D);

	return 0;
}


//svd 
const int MAX_ITER = 100000;
const float eps = 1e-7;

//求取数组x的L2范数
float get_norm(float *x, int n)
{
	float r = 0;
	for (int i = 0; i<n; i++)
		r += x[i] * x[i];
	return sqrt(r);
}

//数组中的每个数都除以其L2范数
float normalize(float *x, int n)
{
	float r = get_norm(x, n);
	if (r<eps)
		return 0;
	for (int i = 0; i<n; i++)
		x[i] /= r;
	return r;
}

//求取矩阵a和矩阵b的内积
inline float product(float*a, float *b, int n)
{
	float r = 0;
	for (int i = 0; i<n; i++)
		r += a[i] * b[i];
	return r;
}

//矩阵规范化
void orth(float *a, float *b, int n)
{
	//|a|=1
	float r = product(a, b, n);
	for (int i = 0; i<n; i++)
		b[i] -= r*a[i];
}

//svd分解矩阵
bool svd(std::vector<std::vector<float>> A, int K, std::vector<std::vector<float>> &U, 
	std::vector<float> &S, std::vector<std::vector<float> > &V)
{
	int M = A.size();
	int N = A[0].size();
	U.clear();
	V.clear();
	S.clear();
	S.resize(K, 0);
	U.resize(K);
	for (int i = 0; i<K; i++)
		U[i].resize(M, 0);

	V.resize(K);
	for (int i = 0; i<K; i++)
		V[i].resize(N, 0);

	//srand(time(0));
	float *left_vector = new float[M];
	float *next_left_vector = new float[M];
	float *right_vector = new float[N];
	float *next_right_vector = new float[N];

	while (1)
	{
		for (int i = 0; i<M; i++)
			left_vector[i] = (float)rand() / RAND_MAX;
		if (normalize(left_vector, M)>eps)
			break;
	}		
	int col = 0;
	for (col = 0; col < K; col++)
	{
		float diff = 1;
		float r = -1;

		for (int iter = 0; diff >= eps && iter < MAX_ITER; iter++)
		{
			memset(next_left_vector, 0, sizeof(float)*M);
			memset(next_right_vector, 0, sizeof(float)*N);
			for (int i = 0; i < M; i++)
				for (int j = 0; j < N; j++)
					next_right_vector[j] += left_vector[i] * A[i][j];
			r = normalize(next_right_vector, N);

			if (r < eps)
				break;

			for (int i = 0; i < col; i++)
				orth(&V[i][0], next_right_vector, N);

			normalize(next_right_vector, N);

			for (int i = 0; i < M; i++)

				for (int j = 0; j < N; j++)

					next_left_vector[i] += next_right_vector[j] * A[i][j];

			r = normalize(next_left_vector, M);

			if (r < eps) 
				break;

			for (int i = 0; i < col; i++)
				orth(&U[i][0], next_left_vector, M);

			normalize(next_left_vector, M);
			diff = 0;
			for (int i = 0; i < M; i++)
			{
				float d = next_left_vector[i] - left_vector[i];
				diff += d*d;
			}

			memcpy(left_vector, next_left_vector, sizeof(float)*M);
			memcpy(right_vector, next_right_vector, sizeof(float)*N);

		}

		if (r >= eps)
		{
			S[col] = r;
			memcpy((char *)&U[col][0], left_vector, sizeof(float)*M);
			memcpy((char *)&V[col][0], right_vector, sizeof(float)*N);
		}
		else
			break;
	}		

	delete[] next_left_vector;
	delete[] next_right_vector;
	delete[] left_vector;
	delete[] right_vector;
	return true;

}

void print(std::vector<std::vector<float>> &A) 
{
	for (int i = 0; i<A.size(); i++) 
	{
		for (int j = 0; j<A[i].size(); j++) 
		{
			std::cout << std::setprecision(3) << A[i][j] << ' ';
		}
		std::cout << std::endl;
	}
}

//求取矩阵的行列式
float MatDet(float *p, int n)
{
	int r, c, m;
	int lop = 0;
	float result = 0;
	float mid = 1;
	if (n != 1)
	{
		lop = (n == 2) ? 1 : n;           
		for (m = 0; m < lop; m++)
		{
			mid = 1;           
			for (r = 0, c = m; r < n; r++, c++)
			{
				mid = mid * (*(p + r*n + c%n));
			}
			result += mid;
		}
		for (m = 0; m < lop; m++)
		{
			mid = 1;            
			for (r = 0, c = n - 1 - m + n; r < n; r++, c--)
			{
				mid = mid * (*(p + r*n + c%n));
			}
			result -= mid;
		}
	}
	else
		result = *p;
	return result;
}


float Creat_M(float *p, int m, int n, int k)
{
	int len;
	int i, j;
	float mid_result = 0;
	int sign = 1;
	float *p_creat, *p_mid;
	len = (k - 1)*(k - 1);            
	p_creat = (float*)calloc(len, sizeof(float)); 
	p_mid = p_creat;
	for (i = 0; i < k; i++)
	{
		for (j = 0; j < k; j++)
		{
			if (i != m && j != n) 
			{
				*p_mid++ = *(p + i*k + j);
			}
		}
	}
	sign = (m + n) % 2 == 0 ? 1 : -1;   
	mid_result = (float)sign*MatDet(p_creat, k - 1);
	free(p_creat);
	return mid_result;
}

bool Gauss(std::vector<std::vector<float>>A, std::vector<std::vector<float>> &B)
{
	int i, j, k;
	float max, temp;
	std::vector<std::vector<float>> t(A);   
	int n = A.size();

	
	B.resize(n);
	for (i = 0; i < n; i++)
	{
		B[i].resize(n,0);
		B[i][i] = 1;
	}	
	
	for (i = 0; i < n; i++)
	{
	
		max = t[i][i];
		k = i;
		for (j = i + 1; j < n; j++)
		{
			if (fabs(t[j][i]) > fabs(max))
			{
				max = t[j][i];
				k = j;
			}
		}

		if (k != i)
		{
			for (j = 0; j < n; j++)
			{
				temp = t[i][j];
				t[i][j] = t[k][j];
				t[k][j] = temp;

				temp = B[i][j];
				B[i][j] = B[k][j];
				B[k][j] = temp;
			}
		}

		if (t[i][i] == 0)
		{
			//cout << "There is no inverse matrix!";
			return false;
		}

		temp = t[i][i];
		for (j = 0; j < n; j++)
		{
			t[i][j] = t[i][j] / temp;        
			B[i][j] = B[i][j] / temp;        
		}
		for (j = 0; j < n; j++)        
		{
			if (j != i)                
			{
				temp = t[j][i];
				for (k = 0; k < n; k++)        
				{
					t[j][k] = t[j][k] - t[i][k] * temp;
					B[j][k] = B[j][k] - B[i][k] * temp;
				}
			}
		}
	}
	//getchar();
	return true;
}

bool CalMatDot(std::vector<std::vector<float>> A, std::vector<std::vector<float>> B, std::vector<std::vector<float>>&C)
{
	int i, j, k;
	int m = A.size();
	int n = A[0].size();
	int mm = B.size();
	int nn = B[0].size();
	
	if (n != mm)
	{
		return false;
	}

	C.clear();
	C.resize(m);
	for (i = 0; i < m; i++)
	{
		C[i].resize(nn, 0);
	}
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < nn; j++)
		{
			for (k = 0; k < n; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return true;
}

/************************************/
bool center_and_normalize_points(std::vector<std::vector<float>> points, std::vector<std::vector<float>> &points_new, 
	std::vector<std::vector<float>> &point_matrix)
{
	int m = points.size();
	int n = points[0].size();
	point_matrix.clear();
	points_new.clear();

	float *centroid = new float[n];
	memset(centroid, 0, sizeof(float)*n);
	int i, j, k;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			centroid[i] += points[j][i];
		}
		centroid[i] /= m;
	}

	float rms = 0.0f;
	float buff = 0.0f;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			buff = points[j][i] - centroid[i];
			rms += buff*buff;
		}		
	}
	rms /= m;
	rms = sqrt(rms);
	
	float norm_factor = sqrt(2) / rms;
	
	point_matrix.resize(3);
	for (i = 0; i < 3; i++)
	{
		point_matrix[i].resize(3, 0);
	}
	point_matrix[0][0] = norm_factor;
	point_matrix[0][2] = -norm_factor * centroid[0];
	point_matrix[1][1] = norm_factor;
	point_matrix[1][2] = -norm_factor * centroid[1];
	point_matrix[2][2] = 1.0f;

	std::vector<std::vector<float>> pointsh(points);
	for (i = 0; i < m; i++)
	{
		pointsh[i].resize(n + 1);
		pointsh[i][n] = 1.0f;		
	}

	std::vector<std::vector<float>>new_pointsh(3, std::vector<float>(m, 0));
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < n+1; j++)
		{
			for (k = 0; k < m; k++)
			{
				new_pointsh[i][k] += point_matrix[i][j] * pointsh[k][j];
			}
		}
	}

	points_new.assign(points.begin(), points.end());
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			points_new[j][i] = new_pointsh[i][j] / new_pointsh[2][j];
		}
	}

	return true;
}


