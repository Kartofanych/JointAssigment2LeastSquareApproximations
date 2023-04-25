//Karim Nasybullin
//k.nasybullin@innopolis.university
#include <bits/stdc++.h>

using namespace std;
#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif
class Matrix{
public:
    float** elements;
    int n = 0;
    int m = 0;

    Matrix(int n){
        this->n = n;
        this->m = n;
        elements = new float*[n];
        for (int i = 0; i < n; i++) {
            elements[i] = new float[m];
        }
    }

    Matrix(int n, int m){
        this->n = n;
        this->m = m;
        elements = new float*[n];
        for (int i = 0; i < n; i++) {
            elements[i] = new float[m];
        }
    }

    void outputMatrix(){
        for(int i = 0; i < n; i ++){
            for(int j = 0; j < m; j ++){
                if(abs(elements[i][j]) < 0.00001){
                    cout<<fixed<<setprecision(4)<<0.0000<<" ";
                }else {
                    cout << fixed << setprecision(4) << elements[i][j] << " ";
                }
            }
            cout<<endl;
        }
    }


    friend istream& operator>>(istream& stream, Matrix& matrix)
    {

        for(int i = 0; i < matrix.n; i ++){
            for(int j = 0; j < matrix.m; j ++){
                stream>>matrix.elements[i][j];
            }
        }

        return stream;
    }

    friend ostream& operator<<(ostream& stream, Matrix& matrix)
    {

        for(int i = 0; i < matrix.n; i ++){
            for(int j = 0; j < matrix.m; j ++){
                stream<<setprecision(2)<<matrix.elements[i][j]<<" ";
            }
            stream<<endl;
        }

        return stream;
    }


    Matrix operator+(Matrix m2){
        auto m3 = Matrix(n, m);

        for(int i = 0; i < m3.n; i ++){
            for(int j = 0; j < m3.m; j ++){
                m3.elements[i][j] = elements[i][j] + m2.elements[i][j];
            }
        }
        return m3;
    }

    Matrix operator-(Matrix m2){
        Matrix m3 = Matrix(n, m);

        for(int i = 0; i < m3.n; i ++){
            for(int j = 0; j < m3.m; j ++){
                m3.elements[i][j] =  elements[i][j] - m2.elements[i][j];
            }
        }
        return m3;
    }

    Matrix operator*(Matrix m2){
        auto m3 = Matrix(n, m2.m);

        for(int i = 0; i < m3.n; i ++){
            for(int j = 0; j < m3.m; j ++){
                float sum = 0;
                for(int k = 0; k < m; k ++){
                    sum += elements[i][k] * m2.elements[k][j];
                }
                m3.elements[i][j] = sum;
            }
        }
        return m3;
    }

    Matrix T(){
        auto m1 = Matrix(m, n);
        for(int i = 0; i < m1.n; i ++){
            for(int j = 0; j < m1.m; j ++){
                m1.elements[i][j] = elements[j][i];
            }
        }
        return m1;
    }


};

class IdentityMatrix : public Matrix{
public:
    IdentityMatrix(int n) : Matrix(n) {
        for(int i = 0; i < n; i ++){
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    elements[i][j] = 1;
                }else{
                    elements[i][j] = 0;
                }
            }
        }
    }
};

class EliminationMatrix : public IdentityMatrix{
private:
    int a = 0;
    int b = 0;
public:
    EliminationMatrix(int n, int a, int b) : IdentityMatrix(n){
        this->a = a;
        this->b = b;
    }
    Matrix ApplyEliminationMatrixFor(Matrix matrix){
        elements[a][b] = -matrix.elements[a][b]/matrix.elements[b][b];
        auto returnMatrix = Matrix(n);
        for(int i = 0; i < n; i ++){
            for(int j = 0; j < n; j ++){
                returnMatrix.elements[i][j] = elements[i][j];
            }
        }
        return returnMatrix;
    }
};

class PermutationMatrix : public IdentityMatrix{
private:
    int a = 0;
    int b = 0;
public:
    PermutationMatrix(int n, int a, int b) : IdentityMatrix(n){
        this->a = a;
        this->b = b;
    }
    Matrix ApplyPermutationMatrixFor(Matrix matrix){
        elements[a][a] = 0;
        elements[b][b] = 0;
        elements[a][b] = 1;
        elements[b][a] = 1;
        auto returnMatrix = Matrix(n);
        for(int i = 0; i < n; i ++){
            for(int j = 0; j < n; j ++){
                returnMatrix.elements[i][j] = elements[i][j];
            }
        }
        return returnMatrix;
    }
};


int maxRow(Matrix matrix, int column_id, int rows, int next_row_id)
{
    float maxm = 0;
    int index = -1;
    for (int row = next_row_id; row < rows; row++)
    {
        if (abs(matrix.elements[row][column_id]) > maxm) {
            index = row;
            maxm = abs(matrix.elements[row][column_id]);
        }
    }
    return index;
}

void outputAugmented(Matrix A, Matrix b){
    for(int i = 0; i < A.n; i ++){
        for(int j = 0; j < A.m; j ++){
            if(abs(A.elements[i][j]) < 0.01){
                cout<<fixed<<setprecision(2)<<0.00<<" ";
            }else {
                cout << fixed << setprecision(2) << A.elements[i][j] << " ";
            }
        }
        for(int j = 0; j < b.m; j ++){
            if(abs(b.elements[i][j]) < 0.01){
                cout<<fixed<<setprecision(2)<<0.00<<" ";
            }else {
                cout << fixed << setprecision(2) << b.elements[i][j] << " ";
            }
        }
        cout<<endl;
    }
}

Matrix inverseMatrix(Matrix a){
    Matrix A = a;
    int n = a.n;
    auto b = Matrix(n, n);
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < n; j ++){
            if(i!=j) {
                b.elements[i][j] = 0;
            }else{
                b.elements[i][j] = 1;
            }
        }
    }


    int step = 1;

    int next_row_id = 0;
    for (int col = 0; col < n; col++)
    {
        int nonzero_row_id = maxRow(A, col, n, next_row_id);
        if (nonzero_row_id >= 0)
        {
            if (nonzero_row_id != next_row_id)
            {
                auto perm = PermutationMatrix(n, next_row_id, nonzero_row_id);
                auto multMatrix = perm.ApplyPermutationMatrixFor(A);
                A = multMatrix*A;
                b = multMatrix*b;
                step++;

                nonzero_row_id = next_row_id;
            }
            for (int row = next_row_id; row < n; row++)
            {
                if (A.elements[row][col] == 0)
                    continue;
                if (row == nonzero_row_id)
                    continue;
                auto elimination = EliminationMatrix(n, row, nonzero_row_id);
                auto multMatrix = elimination.ApplyEliminationMatrixFor(A);
                A = multMatrix*A;
                b = multMatrix*b;
                step++;
            }
            next_row_id++;
        }
        else
        {
            continue;
        }
    }

    for (int col = n-1; col >= 0; col--){
        for (int row = n - 1; row >= 0; row--){
            if (A.elements[row][col] == 0)
                continue;
            if (row >= col)
                continue;
            auto elimination = EliminationMatrix(n, row, col);
            auto multMatrix = elimination.ApplyEliminationMatrixFor(A);
            A = multMatrix*A;
            b = multMatrix*b;
            step++;
        }
    }


    for(int i = 0; i < n; i ++){
        float f = A.elements[i][i];
        A.elements[i][i] = 1;
        for(int j = 0; j < n; j ++){
            b.elements[i][j]/=f;
        }
    }
    return b;
}

int main() {

    int n, st;
    cin>>n;
    int t[n];
    int b[n];
    auto rand = Matrix(n,1);
    for(int i = 0; i < n; i ++){
        cin>>t[i]>>b[i];
        rand.elements[i][0] = ::rand() % 4 + i;
    }
    cin>>st;

    auto A = Matrix(n, st+1);
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < st+1; j ++){
            if(j == 0){
                A.elements[i][j] = 1;
            }else{
                A.elements[i][j] = pow(t[i], j);
            }
        }
    }




    auto B = Matrix(n, 1);
    for(int i = 0; i < n; i ++){
        B.elements[i][0] = b[i];
    }
    cout<<"A:"<<endl;
    A.outputMatrix();

    cout<<"A_T*A:"<<endl;
    (A.T()*A).outputMatrix();

    cout<<"(A_T*A)^-1:"<<endl;
    inverseMatrix((A.T()*A)).outputMatrix();

    cout<<"A_T*b:"<<endl;
    (A.T()*B).outputMatrix();

    cout<<"x~:"<<endl;
    Matrix ans = ((inverseMatrix((A.T()*A)))*(A.T()*B));
    ans.outputMatrix();





    //https://github.com/AlekseyKorshuk/least-square-approximation/
    //helped a lot with gnuplot
    FILE* pipe = _popen(GNUPLOT_NAME, "w");

    string func;
    for (int i = 0; i <= st; i++){
        stringstream ss;
        ss << ans.elements[i][0] ;
        string s;
        ss >> s;
        func+= s;
        func+= "*x**";
        stringstream ss2;
        ss2 << i ;
        string s2;
        ss2 >> s2;
        func+= s2;
        func+=" + ";
    }
    func+="0";
    cout << func << endl;

    int m = func.length();

    char char_array[n + 1];

    strcpy(char_array, func.c_str());

    fprintf(pipe, "set yrange [0:%d]\n", m + 6);
    fprintf(pipe, "f(x) = %s\n", char_array);

    fprintf(pipe, "plot [0:%d] f(x) title 'appr', '-' using 1:2 title 'exp' with points pointtype 5 pointsize 1\n", m);
    for (int i = 0; i < m; i++){
        fprintf(pipe, "%i %f\n", t[i], rand.elements[i][0]);
    }

    fprintf(pipe, "%s\n", "e");

    _pclose(pipe);

}