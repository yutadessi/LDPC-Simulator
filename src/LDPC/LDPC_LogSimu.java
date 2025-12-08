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
        String fileNAMEME = "New16-(6-10)";
        //------------------------------
        String fileNames = fileNAMEME + "-result.csv"; //結果保存ファイル名
        String filePath = fileNAMEME + "-HMatrix.txt"; //検査行列保存ファイル名

        //符号パラメータ
        int n = 1024; //符号長
        int wr = 16; //行重み(n % wr = 0)
        int[] wc = {6,7,8,9,10}; //列重み
        int maxL = 50; //最大反復回数
        int numFrames = 10000; //フレーム数

        //通信路誤り率eの集合
        double[] eValues = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
//            double[] e = {0.05};

        //出力用保存配列-各誤り率のデータ
        double[][][] channelBitErrorRate = new double[wc.length][eValues.length][numFrames]; //各フレームの実際の通信路誤り率
        double[][] sumChannelBitError = new double[wc.length][eValues.length]; //各誤り率における実際の通信路誤り率の合計
        double[][] aveChannelBitErrorRate = new double[wc.length][eValues.length]; //各誤り率における実際の通信路誤り率の平均
        double[][] varianceChannelBitError = new double[wc.length][eValues.length]; //各誤り率における実際の通信路誤り率の分散
        double[][] frameErrorRate = new double[wc.length][eValues.length]; //FER(Frame Error Rate)
        double[][] infoBitErrorRate = new double[wc.length][eValues.length]; //IBER(Info Bit Error Rate)
        double[][] averageTrueIterations = new double[wc.length][eValues.length]; //訂正成功時の平均繰り返し回数
        double[][] averageFalseIterations = new double[wc.length][eValues.length]; //訂正失敗時の平均繰り返し回数
        int[][] residualsErrorBits = new int[wc.length][eValues.length]; //情報ビットの残留誤りビット数
        int[][] errorCorrectionBits = new int[wc.length][eValues.length]; //情報ビットの誤訂正ビット数
        double[][][] iterationDistribution = new double[wc.length][eValues.length][maxL]; //反復回数の度数分布

        //各列重みでのシミュレーション実行
        for(int column = 0;column < wc.length;column++){

            //検査行列Hと生成行列Gの作成
            int [][] h = GenerateMatrix.gallagerCheckMatrix(n,wr,wc[column]);
            int [][] g = GenerateMatrix.generatorMatrix(h);

            //検査行列を保存
            CheckMatrixIO.saveCheckMatrix(h,filePath);

            //HとGを組織符号化
            List<Integer> columnIndicatesToSwap = new ArrayList<>();  //組織符号化用インデックス
            int[][] encodedG = TissueEncoder.EncodeG(g,columnIndicatesToSwap);
            int[][] encodedH = TissueEncoder.EncodeH(h,columnIndicatesToSwap);

            //各通信路誤り率でのシミュレーション
            for(int errorRate = 0; errorRate < eValues.length; errorRate++){

                //正誤毎の反復回数,フレーム数の合計([0]は反復回数,[1]はフレーム数)
                int[] trueInfo = new int[2];
                int[] falseInfo = new int[2];

                for(int frame = 0;frame < numFrames;frame++){
                     //メッセージと送信語、受信語の作成
                    int[] c = GenerateC.geneC(encodedG);
                    int[] r = Channel.GenerateR(c,eValues[errorRate]);

                    //実際の通信路での誤り率の取得
                    channelBitErrorRate[column][errorRate][frame] = Channel.CheckError(c,r);

                    //情報ビットの非誤りビットのインデックス
                    List<Integer> noErrorBitIndex = new ArrayList<>();
                    for(int i = 0;i < g.length;i++){
                        if(c[i] == r[i])noErrorBitIndex.add(i);
                    }

                    //対数領域sum-product復号,確率領域sum-product復号法,Min-Sum復号法
                    LogDecoder.DecodeResult result = LogDecoder.decode(encodedH,r,eValues[errorRate],maxL);
//                    ProbDecoder.DecodingResult result = ProbDecoder.decode(encodedH,r,eValues[errorRate],maxL);
//                    MinSumDecoder.DecodeResult result = MinSumDecoder.decode(encodedH,r,eValues[errorRate],maxL);

                    //復号後と反復回数の取得
                    int[] decodedC = result.decodedCode();
                    int iterations = result.iterationNum();

                    //フレームの正誤判定
                    boolean isFrameCorrect = Arrays.equals(c,decodedC); //trueなら成功.falseなら失敗

                    //反復回数の保存
                    iterationDistribution[column][errorRate][iterations-1] ++;

                    //正誤毎の反復回数,フレーム数の加算
                    if(isFrameCorrect){
                        trueInfo[0] += iterations;
                        trueInfo[1] ++;
                    }else{
                        falseInfo[0] += iterations;
                        falseInfo[1] ++;
                    }

                    //残留誤りビットと誤訂正ビットの加算
                    for(int z = 0;z < g.length;z++){
                        if(c[z] != decodedC[z])residualsErrorBits[column][errorRate] ++;
                    }
                    for(int x : noErrorBitIndex){
                        if(c[x] != decodedC[x])errorCorrectionBits[column][errorRate] ++;
                    }

                    //実際の誤り率の加算
                    sumChannelBitError[column][errorRate] += channelBitErrorRate[column][errorRate][frame];
                }
                //正誤毎の反復回数の平均
                averageTrueIterations[column][errorRate] = (double)trueInfo[0]/ trueInfo[1];
                averageFalseIterations[column][errorRate] = (double)falseInfo[0]/ falseInfo[1];

                //FER
                frameErrorRate[column][errorRate] = (double)falseInfo[1]/numFrames;

                //IBER
                infoBitErrorRate[column][errorRate] = (double)residualsErrorBits[column][errorRate] / (g.length * numFrames);

                //実際の誤り率の平均
                aveChannelBitErrorRate[column][errorRate] = sumChannelBitError[column][errorRate]/numFrames;
            }
        }

        //実際の通信路誤り率の分散を求める
        for(int i = 0;i < wc.length;i++){
            for(int j = 0;j < eValues.length;j++){
                for(int k = 0;k < numFrames;k++){
                    varianceChannelBitError[i][j] += Math.pow((channelBitErrorRate[i][j][k] - aveChannelBitErrorRate[i][j]),2);
                }
                varianceChannelBitError[i][j] /= numFrames;
            }
        }

        try (PrintWriter pw = new PrintWriter(fileNames, Charset.forName("Windows-31j"))){
            for(int i = 0;i < wc.length;i++)pw.printf("%s,%s,%s,%s,%s,%s,,,,","符号長","行重み","列重み","符号化率","最大反復回数","フレーム数");
            pw.printf("\n");
            for(int i = 0;i < wc.length;i++) pw.printf("%s,%s,%s,%s,%s,%s,,,,",n,wr,wc[i],(1-(double)wc[i]/wr),maxL,numFrames);
            pw.printf("\n\n");
            for(int i = 0;i < wc.length;i++)pw.printf("%s,%s,%s,%s,%s,%s,%s,%s,,","通信路誤り率","実際の通信路誤り率の平均","実際の通信路誤り率の分散","FER","IBER","成功時の平均繰り返し回数","失敗時の平均繰り返し回数","平均誤訂正ビット率");
            pw.printf("\n");
            for(int i = 0;i < eValues.length;i++){
                for(int j = 0;j < wc.length;j++){
                    pw.printf("%.2f,%s,%s,%.4f,%s,%s,%s,%s,,", eValues[i], aveChannelBitErrorRate[j][i],varianceChannelBitError[j][i], frameErrorRate[j][i], (residualsErrorBits[j][i]/numFrames),averageTrueIterations[j][i],averageFalseIterations[j][i], (errorCorrectionBits[j][i]/residualsErrorBits[j][i]));
                }
                pw.printf("\n");
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
