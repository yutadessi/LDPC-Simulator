package LDPC;

import java.util.ArrayList;
import java.util.List;

public class LogDecoder {

    //逆ハイパボリックタンジェント
    public static double atanh(double x){
        return 0.5 * Math.log( ( 1 + x ) / ( 1 - x ) );
    }

    //レコード
    public record DecodeResult (int[] decodedCode, int iterationCount){}

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

        //基本処理
        for(int l = 0;l < maxL;l++){

            //変数ノード処理
            for(int j = 0;j < numV;j++){
                List<Integer> connectedC = B.get(j);
                for(int i : connectedC){

                }
            }
        }


    }

}
