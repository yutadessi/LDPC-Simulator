package LDPC;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;

public class LoadH_Prob {
    public static void main(String[] args) {

        //ファイル名、毎回変える！！--------
        String fileNAMEME = "No.3";
        //------------------------------

        String fileNames = fileNAMEME + "-LoadHResult.txt";
        String filePath = fileNAMEME + "-HMatrix.txt";
        try (PrintWriter pw = new PrintWriter(fileNames, StandardCharsets.UTF_8)){

            //符号パラメーター
            int maxL = 20; //最大反復回数

            //シミュレーション設定
            int numFrames = 1000;

            //通信路誤り率eの設定
//            double [] eValues = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
            double[] eValues = {0.03};

            //検査行列Hと生成行列Gの作成
            int [][] H = CheckMatrixIO.loadCheckMatrix(filePath);
            int [][] G = GenerateMatrix.generatorMatrix(H);

            int n = H[0].length;
            int wr = 0;
            for(int i = 0;i < n;i++){
                if(H[0][i] == 0){
                    wr = i;
                    break;
                }
            }
            int wc = (H.length * wr / n);

            pw.println("n=" + n + ",wr=" + wr + ",wc=" + wc + ",maxL=" + maxL);

            //HとGを組織符号化
            List<Integer> columnIndicatesToSwap = new ArrayList<>();
            int[][] encodedG = TissueEncoder.EncodeG(G,columnIndicatesToSwap);
            int[][] encodedH = TissueEncoder.EncodeH(H,columnIndicatesToSwap);

            //フレーム毎の表示達
            double[][] groupOfCBER = new double[eValues.length][numFrames];
            String[][] groupOfFrame = new String[eValues.length][numFrames];
            long[][] groupOfErroIBits = new long[eValues.length][numFrames];
            int[][] groupOfIterations = new int[eValues.length][numFrames];

            pw.printf("%8s | %13s | %12s\n","E","FER","IBER");

            int num = 0;

            //復号
            for(double e : eValues){
                long frameErrorCount = 0;
                long frameErrorCountPerFrame = 0;
                long bitErrorCount = 0;
                long totalInfoBits = 0;
                int infoBitLength = encodedG.length;
                long BECountPerFrame = 0;

                for(int frame = 0;frame < numFrames;frame++){

                    //メッセージと送信語の作成
                    int[] c = GenerateC.geneC(encodedG);

                    //受信語作成
                    int[] r = Channel.GenerateR(c,e);

                    double cBER = Channel.CheckError(c,r);

                    //確率領域sum-product復号
                    ProbDecoder.DecodingResult result = ProbDecoder.decode(encodedH,r,e,maxL);

                    int[] estimatedC = result.decodedCode();
                    int iterations = result.iterationNum();

                    //フレーム誤りカウント
                    if(!Arrays.equals(c,estimatedC)){
                        frameErrorCount++;
                    }
                    String frameError = (frameErrorCount-frameErrorCountPerFrame) > 0 ? "False" : "True";

                    //情報ビット誤りカウント
                    totalInfoBits += infoBitLength;

                    for(int i = 0;i < infoBitLength;i++){
                        if(c[i] != estimatedC[i]) bitErrorCount++;
                    }
                    groupOfCBER[num][frame] = cBER;
                    groupOfFrame[num][frame] = frameError;
                    groupOfErroIBits[num][frame] = (bitErrorCount - BECountPerFrame);
                    groupOfIterations[num][frame] = iterations;

                    frameErrorCountPerFrame = frameErrorCount;
                    BECountPerFrame = bitErrorCount;
                }
                double fer = (double)frameErrorCount/numFrames;
                double iber = (double)bitErrorCount/totalInfoBits;
                pw.printf("%6.2f | %10.6f | %12.8f\n", e, fer, iber);
                num++;
            }

            //各フレームの情報表示
            pw.printf("%-13s | %-5s | %-10s | %-6s\n","C-BER","Frame","ErrorIBits","Iterations");
            for(int i = 0;i < eValues.length;i++){
                for(int j = 0;j < numFrames;j++){
                    pw.printf("%10.8f | %6s | %12d | %6d\n",groupOfCBER[i][j],groupOfFrame[i][j],groupOfErroIBits[i][j],groupOfIterations[i][j]);
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
