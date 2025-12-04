package LDPC;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.IOException;
import java.nio.charset.Charset;

//--------------------デスクトップPCでの処理速度--------------------
//処理速度(フレーム数:10000,1024-8-4サイズ,誤り率数:10,Lmax:20 ):18分
//処理速度(フレーム数:10000,1024-8-4サイズ,誤り率数:10,Lmax:100):31分

public class LDPC_LogSimu {
    public static void main(String[] args) {

        //ファイル名、毎回変える！！--------
        String fileNAMEME = "Trial";
        //------------------------------

        String fileNames = fileNAMEME + "-result.csv";
        String filePath = fileNAMEME + "-HMatrix.txt";
        try (PrintWriter pw = new PrintWriter(fileNames, Charset.forName("Windows-31j"))){

            //符号パラメーター
            int n = 1024; //符号長
            int wr = 8; //行重み,n % wr == 0
            int wc = 3; //列重み
            int maxL = 50; //最大反復回数

            //シミュレーション設定
            int numFrames = 10;

            pw.printf("%s,%s,%s,%s,%s\n",n,wr,wc,maxL,numFrames);

            //通信路誤り率eの設定
            double[] eValues = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
//            double[] eValues = {0.05};

            //検査行列Hと生成行列Gの作成
            int [][] H = GenerateMatrix.gallagerCheckMatrix(n,wr,wc);
            int [][] G = GenerateMatrix.generatorMatrix(H);

            //検査行列を保存
            CheckMatrixIO.saveCheckMatrix(H,filePath);

            //HとGを組織符号化
            List<Integer> columnIndicatesToSwap = new ArrayList<>();  //組織符号化用インデックス
            int[][] encodedG = TissueEncoder.EncodeG(G,columnIndicatesToSwap);
            int[][] encodedH = TissueEncoder.EncodeH(H,columnIndicatesToSwap);

            //フレーム毎の表示達
            double[][] groupOfCBER = new double[eValues.length][numFrames];
            String[][] groupOfFrame = new String[eValues.length][numFrames];
            long[][] groupOfErrorInfoBits = new long[eValues.length][numFrames];
            int[][] groupOfIterations = new int[eValues.length][numFrames];

            pw.printf("%s,%s,%s,%s,%s,%s\n","通信路誤り率","実際の通信路誤り率の平均","FER","IBER","成功時の平均繰り返し回数","失敗時の平均繰り返し回数");

            int num = 0;

            //復号
            for(double e : eValues){
                long frameErrorCount = 0; //フレーム誤り数合計
                long frameErrorCountPreviousTotal = 0; //前フレームまでのフレーム誤り数合計
                long bitErrorCount = 0; //ビットエラー
                long totalInfoBits = 0;
                int infoBitLength = encodedG.length;
                long BECountPerFrame = 0;

                double aveTrueIterations;
                double aveFalseIterations;
                int[] sumTrueIterations = new int[2];
                int[] sumFalseIterations = new int[2];

                //各誤り率での実際の通信路誤り率の合計と平均
                double aveCBER;
                double sumCBER = 0;

                for(int frame = 0;frame < numFrames;frame++){

                    //メッセージと送信語の作成
                    int[] c = GenerateC.geneC(encodedG);

                    //受信語作成
                    int[] r = Channel.GenerateR(c,e);

                    double cBER = Channel.CheckError(c,r);

                    //対数領域sum-product復号
                    LogDecoder.DecodeResult result = LogDecoder.decode(encodedH,r,e,maxL);

                    int[] estimatedC = result.decodedCode();
                    int iterations = result.iterationNum();

                    //フレーム誤りカウント
                    if(!Arrays.equals(c,estimatedC)){
                        frameErrorCount++;
                    }
                    String frameError = (frameErrorCount-frameErrorCountPreviousTotal) > 0 ? "False" : "True";

                    //情報ビット誤りカウント
                    totalInfoBits += infoBitLength;

                    for(int i = 0;i < infoBitLength;i++){
                        if(c[i] != estimatedC[i]) bitErrorCount++;
                    }
                    groupOfCBER[num][frame] = cBER;
                    groupOfFrame[num][frame] = frameError;
                    groupOfErrorInfoBits[num][frame] = (bitErrorCount - BECountPerFrame);
                    groupOfIterations[num][frame] = iterations;

                    frameErrorCountPreviousTotal = frameErrorCount;
                    BECountPerFrame = bitErrorCount;

                    //実際の誤り率の合計
                    sumCBER += cBER;

                    //正誤毎の平均繰り返し回数
                    if(frameError == "True"){
                        sumTrueIterations[0] += iterations;
                        sumTrueIterations[1] ++;
                    }else{
                        sumFalseIterations[0] += iterations;
                        sumFalseIterations[1] ++;
                    }

                }
                double fer = (double)frameErrorCount/numFrames;
                double iber = (double)bitErrorCount/totalInfoBits;

                aveCBER = sumCBER / numFrames;
                aveTrueIterations = (double)sumTrueIterations[0] / sumTrueIterations[1];
                aveFalseIterations = (double)sumFalseIterations[0] / sumFalseIterations[1];

                pw.printf("%.2f,%s,%.4f,%s,%s,%s\n", e, aveCBER, fer, iber,aveTrueIterations,aveFalseIterations);
                num++;
            }

            pw.println(" ");

            //各フレームの情報表示
            pw.printf("%s,%s,%s,%s\n","C-BER","Frame","ErrorIBits","Iterations");
            for(int i = 0;i < eValues.length;i++){
                for(int j = 0;j < numFrames;j++){
                    pw.printf("%s,%s,%s,%s\n",groupOfCBER[i][j],groupOfFrame[i][j],groupOfErrorInfoBits[i][j],groupOfIterations[i][j]);
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
