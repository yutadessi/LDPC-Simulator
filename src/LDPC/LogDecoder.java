package LDPC;

import java.util.ArrayList;
import java.util.List;

public class LogDecoder {

    //逆ハイパボリックタンジェント
    public static double atanh(double x){
        return 0.5 * Math.log( ( 1.0 + x ) / ( 1.0 - x ) );
    }


    //レコード
    public record DecodeResult (int[] decodedCode, int iterationNum){}

    public static DecodeResult decode (int[][] encodedH, int[] r, double e, int maxL){

        //初期条件
        int numC = encodedH.length;
        int numV = encodedH[0].length;

         //インデックスA,B
        List<List<Integer>> A = new ArrayList<>(numC); // A(i)
        for (int i = 0; i < numC; i++) {
            A.add(new ArrayList<>());
        }
        List<List<Integer>> B = new ArrayList<>(numV); // B(j)
        for (int j = 0; j < numV; j++) {
            B.add(new ArrayList<>());
        }
        for (int i = 0; i < numC; i++) {
            for (int j = 0; j < numV; j++) {
                if (encodedH[i][j] == 1) {
                    B.get(j).add(i); // B(j) に i を追加
                    A.get(i).add(j); // A(i) に j を追加
                }
            }
        }

         //対数尤度比λ
        double[] lamda = new double[r.length];
        double l0 = Math.log((1 - e) / e);

        for(int j = 0;j < r.length;j++){
            lamda[j] = (r[j] == 0) ? l0 : -l0;
        }

         //対数領域メッセージα,β
        double[][] alpha = new double[numC][numV];
        double[][] beta = new double[numV][numC];

        //推定語
        int[] estimatedCode = new int[numV];

        //基本処理
        for(int l = 0;l < maxL;l++){

            //変数ノード処理
            for(int j = 0;j < numV;j++){
                List<Integer> connectedC = B.get(j);
                for(int i : connectedC){
                    double product = lamda[j];//λ(j)

                    for(int k : connectedC){//+ Σα
                        if(k == i) continue;
                        product += alpha[k][j];
                    }

                    beta[j][i] = product;//代入
                }
            }

            //チェックノード処理
            for(int i = 0;i < numC;i++){
                List<Integer> connectedV = A.get(i);

                for(int j : connectedV){
                    double product = 1.0;

                    for(int k : connectedV){//Πtanh(0.5*β)
                        if(k == j) continue;
                        product *= Math.tanh(0.5 * beta[k][i]);
                    }

                    alpha[i][j] = 2 * atanh(product);// 2tanh^-1
                }
            }

            //一時推定ビットの決定
            double[] gamma = new double[numV];
            for(int j = 0;j < numV;j++){
                gamma[j] = lamda[j];
                for(int k : B.get(j)){
                    gamma[j] += alpha[k][j];
                }
                estimatedCode[j] = (gamma[j] < 0)? 1 : 0;
            }

            //パリティ検査
            int sumSyndro = 0;
            int[] syndrome = new int[encodedH.length];
            for(int m = 0;m < encodedH.length;m++){
                for(int n = 0;n < encodedH[0].length;n++){
                    syndrome[m] += encodedH[m][n] * estimatedCode[n];
                }
                syndrome[m] %= 2;
                sumSyndro += syndrome[m];
            }
            if(sumSyndro == 0){
                return new DecodeResult(estimatedCode,l + 1);
            }
        }

        return new DecodeResult(estimatedCode,maxL);
    }

}
