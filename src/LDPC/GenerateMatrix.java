package LDPC;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

//行列作成クラス
public class GenerateMatrix {
    //検査行列作成メソッド
    public static int[][] gallagerCheckMatrix(int n,int wr,int wc){
        int m = n * wc / wr;//検査行列の行数
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
    public static int[][] generatorMatrix(int[][] H){
        int n = H[0].length;
        int m = H.length;
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

        //生成行列Gの抽出
        for(int i = 0;i < n - m;i++){
            for(int j = 0;j < n;j++){
                G_Matrix[i][j] = HtI[m + i][m + j];
            }
        }
        return G_Matrix;
    }

    public static int[][] generateQC(int z,int mb,int nb){
        int[][] baseM = new int[mb][nb];
        Random rand = new Random();
        boolean fails;

        do{
            //ベース行列の作成
            for(int i = 0;i < mb;i++){
                for(int j = 0;j < nb;j++){
                    baseM[i][j] = rand.nextInt(z+1) - 1;
                }
            }

            fails = false;

            // --- 条件1：各行の非-1要素 ≥ 2 ---
            for (int i = 0; i < mb; i++) {
                int count = 0;
                for (int j = 0; j < nb; j++) {
                    if (baseM[i][j] != -1) count++;
                }
                if (count < 2) {
                    fails = true;
                    break;
                }
            }

            // --- 条件2：各列の非-1要素 ≥ 2 ---
            if (!fails) {
                for (int j = 0; j < nb; j++) {
                    int count = 0;
                    for (int i = 0; i < mb; i++) {
                        if (baseM[i][j] != -1) count++;
                    }
                    if (count < 2) {
                        fails = true;
                        break;
                    }
                }
            }

            // --- 条件3：4-cycle チェック ---
            if (!fails && has4Cycle(baseM, z)) {
                fails = true;
            }

        }while (fails);

        //展開
        int[][] QC = new int[mb * z][nb * z];
        for(int i = 0;i < mb;i++){
            for(int j = 0;j < nb;j++){
                if(baseM[i][j] == -1) continue;
                for(int k = 0;k < z;k++){
                    int r = i * z + k;
                    int c = j * z + (k + baseM[i][j]) % z;
                    QC[r][c] = 1;
                }
            }
        }

        return QC;
    }

    static boolean has4Cycle(int[][] B, int z) {
        int mb = B.length;
        int nb = B[0].length;

        for (int i1 = 0; i1 < mb; i1++) {
            for (int i2 = i1 + 1; i2 < mb; i2++) {
                for (int j1 = 0; j1 < nb; j1++) {
                    if (B[i1][j1] == -1 || B[i2][j1] == -1) continue;
                    for (int j2 = j1 + 1; j2 < nb; j2++) {
                        if (B[i1][j2] == -1 || B[i2][j2] == -1) continue;

                        int lhs = (B[i1][j1] - B[i1][j2] + z) % z;
                        int rhs = (B[i2][j1] - B[i2][j2] + z) % z;

                        if (lhs == rhs) {
                            return true; // 4-cycle あり
                        }
                    }
                }
            }
        }
        return false;
    }

}
