package LDPC;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;

public class LDPC_SimuGT {
    public static void main(String[] args) {

        String fileNames = "result-001-01().txt";
        try (PrintWriter pw = new PrintWriter(fileNames, StandardCharsets.UTF_8)){

            //符号パラメーター
            int n = 1024; //符号長
            int wr = 8; //行重み,n % wr == 0
            int wc = 4; //列重み
            int maxL = 20; //最大反復回数

            //シミュレーション設定
            int numFrames = 1;

            pw.println("n=" + n + ",wr=" + wr + ",wc=" + wc + ",maxL=" + maxL);

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

            //復号
            for(double e : eValues){
                long frameErrorCount = 0;
                long bitErrorCount = 0;
                long totalInfoBits = 0;
                int infoBitLength = encodedG.length;

                for(int frame = 0;frame < numFrames;frame++){

                    //メッセージと送信語の作成
                    int[] c = GenerateC.geneC(encodedG);

                    //受信語作成
                    int[] r = Channel.GenerateR(c,e);

                    //確率領域sum-product復号
                    ProbDecoder.DecodingResult result = ProbDecoder.decode(encodedH,r,e,maxL,pw);

                    int[] estimatedC = result.decodedCodeword;
                    int iterations = result.iterationCount;

                    //フレーム誤りカウント
                    if(!Arrays.equals(c,estimatedC)){
                        frameErrorCount++;
                    }

                    //情報ビット誤りカウント
                    totalInfoBits += infoBitLength;

                    for(int i = 0;i < infoBitLength;i++){
                        if(c[i] != estimatedC[i]) bitErrorCount++;
                    }
                }
                double fer = (double)frameErrorCount/numFrames;
                double iber = (double)bitErrorCount/totalInfoBits;
                pw.printf("%.2f \t\t| %.6f \t| %.8f\n", e, fer, iber);
            }
            Print.DetailMatrix(encodedH,pw);//書き写す用
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
