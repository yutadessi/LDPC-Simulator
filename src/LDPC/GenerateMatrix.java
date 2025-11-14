package LDPC;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

//行列作成クラス
public class GenerateMatrix {
    //検査行列作成メソッド
    public static int[][] gallagerCheckMatrix(int n,int wr,int wc){
        int m = n*wc/wr;//検査行列の行数
        int blockRows = n/wr;//部分行列の行数
        int[][] H_matrix = new int [m][n];

        //最初の部分行列作成
        for(int i = 1;i <= n/wr;i++){
            for(int j = 0;j < n;j++){
                H_matrix[(i-1)][j] = j >= (i-1)*wr && j < i*wr ? 1 : 0;
            }
        }
        //以降の部分行列作成
        List<Integer> index = new ArrayList<>();
        for(int i = 0;i < n;i++) index.add(i);
        for(int blockNum = 1;blockNum < wc;blockNum++){//blockNumは現在の部分行列の番号
            Collections.shuffle(index);//列ベクトルのランダム置換
            for(int i = 0;i < n;i++){
                for(int j = 0;j < blockRows;j++){
                    H_matrix[blockRows*blockNum+j][i] = H_matrix[j][index.get(i)];
                }
            }
        }
        return H_matrix;
    }

    //生成行列作成
    public static int[][] generatorMatrix(int[][] H,int n,int wr,int wc){
        int m = n * wc / wr;
        int[][] HtI = new int[n][n + m];//拡大行列HtI
        int[][] G_Matrix = new int[n - m][n];//生成行列G
        //拡大行列作成
        for(int i = 0;i < n;i++){ //検査行列の転置
            for(int j = 0;j < m;j++){
                HtI[i][j] = H[j][i];
            }
        }
        for(int i = 0;i < n;i++){ //単位行列作成
            for(int j = m;j < m + n;j++){
                HtI[i][j] = i + m == j ? 1 : 0;
            }
        }

        //ガウスの消去法
        for(int i = 0;i < m;i++){//i列目に1があるか確認、交換、加算をここで行う。
            List<Integer> tempo = new ArrayList<>();
            int[] temp;
            for (int j = 0;j < n;j++){//1発見器
                if(j > i && HtI[j][i] == 1) tempo.add(j);
            }
            if(tempo.size() == 0)continue;
            if(HtI[i][i] != 1){//交換
                temp = HtI[i];
                HtI[i] = HtI[tempo.get(0)]; 
                HtI[tempo.get(0)] = temp;
            }
            for(int j = 0;j < tempo.size();j++){//加算
                if(HtI[tempo.get(j)][i] == 1){
                    for(int k = 0;k < n + m;k++){
                        if(HtI[i][k] == 1 && HtI[tempo.get(j)][k] == 0){
                            HtI[tempo.get(j)][k] = 1;
                            continue;
                        }
                        if(HtI[i][k] == 1 && HtI[tempo.get(j)][k] == 1)HtI[tempo.get(j)][k] = 0;
                    }
                }
            }
            tempo.clear();
        }
//        Print.Matrix(HtI);

        //生成行列Gの抽出
        for(int i = 0;i < n - m;i++){
            for(int j = 0;j < n;j++){
                G_Matrix[i][j] = HtI[m + i][m + j];
            }
        }
        return G_Matrix;
    }
}
