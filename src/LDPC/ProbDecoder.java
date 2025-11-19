package LDPC;

import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.io.PrintWriter;

public class ProbDecoder {
    public static class DecodingResult {
        public final int[] decodedCodeword;
        public final int iterationCount;

        public DecodingResult(int[] decodedCodeword,int iterationCount){
            this.decodedCodeword = decodedCodeword;
            this.iterationCount = iterationCount;
        }
    }
    public static DecodingResult decode(int[][] encodedH, int[] r, double e, int maxL,PrintWriter pw) {

        //各ノードサイズ
        int numC = encodedH.length;
        int numV = encodedH[0].length;

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

        //尤度
        double[][] likelihoods = new double[numV][2];
        double p_correct = 1.0 - e;
        double p_error = e;

        for (int j = 0; j < numV; j++) {
            if (r[j] == 0) {
                likelihoods[j][0] = p_correct;
                likelihoods[j][1] = p_error;
            } else {
                likelihoods[j][0] = p_error;
                likelihoods[j][1] = p_correct;
            }
        }

        //初期化
        double[][][] msgCtoV = new double[numC][numV][2];
        for (int i = 0; i < numC; i++) {
            for (int j : A.get(i)) {
                msgCtoV[i][j][0] = 1.0;
                msgCtoV[i][j][1] = 1.0;
            }
        }

        int[] estimatedC = new int[numV];

        //基本処理
        for (int l = 0; l < maxL; l++) {
            double[][][] msgVtoC = new double[numV][numC][2];

            //変数ノード処理
            for (int j = 0; j < numV; j++) {
                List<Integer> connectedC = B.get(j);
                for (int i : connectedC) {
                    double p0 = likelihoods[j][0];
                    double p1 = likelihoods[j][1];

                    for (int k : connectedC) {
                        if (k != i) {
                            p0 *= msgCtoV[k][j][0];
                            p1 *= msgCtoV[k][j][1];
                        }
                    }
                    msgVtoC[j][i][0] = p0;
                    msgVtoC[j][i][1] = p1;
                }
            }

            //チェックノード処理
            for(int i = 0;i < numC;i++){
                List<Integer> connectedV = A.get(i);

                for(int j : connectedV){
                    int[] vNodes = new int[connectedV.size() - 1];
                    double sum0 = 0,sum1 = 0;
                    long svmax = 1L << vNodes.length;
                    for(int sigma = 0;sigma < svmax;sigma++){
                        double product = 1.0;

                        //Vi∈Ai\j
                        for(int countL = 0;countL < vNodes.length;countL++){
                            vNodes[countL] = (int)((sigma >> countL) & 1);
                        }

                        //Π
                        for(int k = 0;k < vNodes.length;k++){
                            if(connectedV.get(k) < j){
                                product *= msgVtoC[connectedV.get(k)][i][vNodes[k]];
                            }else if(connectedV.get(k) >= j){
                                product *= msgVtoC[connectedV.get(k+1)][i][vNodes[k]];
                            }
                        }

                        //I
                        if(Arrays.stream(vNodes).sum() % 2 == 0){
                            sum0 += product;
                        }else if(Arrays.stream(vNodes).sum() % 2 == 1){
                            sum1 += product;
                        }
                    }
                    msgCtoV[i][j][0] = sum0;
                    msgCtoV[i][j][1] = sum1;
                }
            }

            //一時推定ビットの決定
            for(int j = 0;j < numV;j++){
                double guess0 = likelihoods[j][0];
                double guess1 = likelihoods[j][1];
                for(int k : B.get(j)){
                    guess0 *= msgCtoV[k][j][0];
                    guess1 *= msgCtoV[k][j][1];
                }
                estimatedC[j] = guess0 >= guess1 ? 0 : 1;
            }

            //パリティ検査
            int sumSyndro = 0;
            int[] syndrome = new int[encodedH.length];
            for(int m = 0;m < encodedH.length;m++){
                for(int n = 0;n < encodedH[0].length;n++){
                    syndrome[m] += encodedH[m][n] * estimatedC[n];
                }
                syndrome[m] %= 2;
                sumSyndro += syndrome[m];
            }
            if(sumSyndro == 0){
//                System.out.println("反復回数" + (l + 1) + "でパリティ検査通過");
                return new DecodingResult(estimatedC,l + 1);
            }

        }
//        System.out.println("最大反復回数" + maxL + "に到達しました。");
        return new DecodingResult(estimatedC,maxL);
    }
}