package LDPC;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;

public class LDPC_SimuGT {
    public static void main(String[] args) {

        String fileNames = "myresult.txt";
        try (PrintWriter pw = new PrintWriter(fileNames, StandardCharsets.UTF_8)){

            //符号パラメーター
            int n = 6; //符号長
            int wr = 2; //行重み
            int wc = 1; //列重み
            int maxL = 2; //最大反復回数

            //シミュレーション設定
            int numFrames = 1;

            //通信路誤り率eの設定
//            double[] eValues = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
            double[] eValues = {0.03};

            //検査行列Hと生成行列Gの作成
            int [][] H = GenerateMatrix.gallagerCheckMatrix(n,wr,wc);
            int [][] G = GenerateMatrix.generatorMatrix(H,n,wr,wc);

            //HとGを組織符号化
            List<Integer> columnIndicatesToSwap = new ArrayList<>();
            int[][] encodedG = TissueEncoder.EncodeG(G,columnIndicatesToSwap);
            int[][] encodedH = TissueEncoder.EncodeH(H,columnIndicatesToSwap);
            pw.println(columnIndicatesToSwap+"\n");
            pw.print("Encoded G ");
            Print.Matrix(encodedG,pw);
            pw.print("Encoded H ");
            Print.Matrix(encodedH,pw);

            //復号
            for(double e : eValues){
                long frameErrorCount = 0;

                for(int frame = 0;frame < numFrames;frame++){

                    //メッセージと送信語の作成と表示
                    int[] c = GenerateC.geneC(encodedG);

                    //受信語作成
                    int[] r = Channel.GenerateR(c,e);

                    //確率領域sum-product復号
                    int[] estimatedC = ProbDecoder.decode(encodedH,r,e,maxL);

                    if(!Arrays.equals(c,estimatedC)){
                        frameErrorCount++;
                    }
                }
                double fer = (double)frameErrorCount/numFrames;
                System.out.printf("%.2f \t\t| %.6f\n",e,fer);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
