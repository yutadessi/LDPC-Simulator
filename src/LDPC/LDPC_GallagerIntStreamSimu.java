package LDPC;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.stream.IntStream;

public class LDPC_GallagerIntStreamSimu {
    public static void main(String[] args) {

        //ファイル名、毎回変える！！--------
        String fileNAMEME = "final8-2(10_000-25)G-Log";
        //------------------------------
        String fileNames = fileNAMEME + "-result.csv";

        //符号パラメータ
        int n = 1024; //符号長
        int wr = 8; //行重み(n % wr = 0)
        int wc = 2; //列重み
        int maxL = 50; //最大反復回数
        int numFrames = 10_000; //フレーム数

        //検査行列作成数
        int numCM = 25;

        //通信路誤り率eの集合
        double[] eValues = {0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05};
//            double[] e = 0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15

        //タイム計測用配列
        double[][] executionTimes = new double[numCM][3]; //トータルの実行時間
        double[][] decodeTimes = new double[numCM][eValues.length]; //各誤り率の復号時間

        //出力用保存配列-各誤り率のデータ
        double[][][] actualChannelBitErrorRate = new double[numCM][eValues.length][numFrames]; //各フレームの実際の通信路誤り率
        double[][] sumChannelBitError = new double[numCM][eValues.length]; //各誤り率における実際の通信路誤り率の合計
        double[][] aveChannelBitErrorRate = new double[numCM][eValues.length]; //各誤り率における実際の通信路誤り率の平均
        double[][] varianceChannelBitError = new double[numCM][eValues.length]; //各誤り率における実際の通信路誤り率の分散
        double[][] frameErrorRate = new double[numCM][eValues.length]; //FER(符号長の失敗確率)
        double[][] informationFrameErrorRate = new double[numCM][eValues.length]; //IFER(情報ビットのみの失敗確率）
        double[][] checkFrameErrorRate = new double[numCM][eValues.length]; //CFER(検査ビットのみの失敗確率）
        double[][] infoBitErrorRate = new double[numCM][eValues.length]; //IBER(Info Bit Error Rate)
        double[][] averageTrueIterations = new double[numCM][eValues.length]; //訂正成功時の平均繰り返し回数
        double[][] averageFalseIterations = new double[numCM][eValues.length]; //訂正失敗時の平均繰り返し回数
        int[][] residualsErrorInfoBits = new int[numCM][eValues.length]; //情報ビットの残留誤りビット数
        double[][] avgErrorCheckBits = new double[numCM][eValues.length]; //検査ビットの残留誤りビット数
        int[][] errorCorrectionBits = new int[numCM][eValues.length]; //情報ビットの誤訂正ビット数
        int[][][] iterationDistribution = new int[numCM][eValues.length][maxL]; //反復回数の度数分布
        int[][] undetectedErrors = new int[numCM][eValues.length]; //シンドロームは0だが,誤訂正している数

        //各列重みでのシミュレーション実行
        for(int column = 0;column < numCM;column++){

            //IntStream用にfinal化
            final int cIndex = column;

            //実行時間全体の計測開始
            long startTotal = System.currentTimeMillis();

            //検査行列hと生成行列gの作成
            int [][] h = GenerateMatrix.gallagerCheckMatrix(n,wr,wc);
            int [][] g = GenerateMatrix.generatorMatrix(h);

            int gLength = g.length;

            //検査行列を保存
            String filePath = fileNAMEME + column + "-HMatrix.txt"; //検査行列保存ファイル名
            CheckMatrixIO.saveCheckMatrix(h,filePath);

            //HとGを組織符号化
            List<Integer> columnIndicatesToSwap = new ArrayList<>();  //組織符号化用インデックス
            int[][] encodedG = TissueEncoder.EncodeG(g,columnIndicatesToSwap);
            int[][] encodedH = TissueEncoder.EncodeH(h,columnIndicatesToSwap);

            //各通信路誤り率でのシミュレーション
            for(int errorRate = 0; errorRate < eValues.length; errorRate++){

                //IntStream用にfinal化
                final int eIndex = errorRate;

                //復号時間の合計用変数を初期化
                final long[] sumDecodeTime = {0};

                //正誤毎の反復回数,フレーム数の合計([0]は反復回数,[1]はフレーム数)
                int[] trueIterations = new int[2];
                int[] falseIterations = new int[2];

                final int[] errorInfoFrameCounter = {0};
                final int[] errorCheckFrameCounter = {0};

                final int[] errorCheckBitsCounter = {0};

                Object lock = new Object();

                //並列処理開始
                IntStream.range(0, numFrames).parallel().forEach(frame -> {

                    //メッセージと送信語、受信語の作成
                    int[] c = GenerateC.geneC(encodedG);
                    int[] r = Channel.GenerateR(c,eValues[eIndex],gLength);

                    //フレームごとの情報ビットの正誤
                    int currentInfoFrameErrorBits = 0;
                    int currentCheckFrameErrorBits = 0;

                    //実際の通信路での誤り率の取得
                    actualChannelBitErrorRate[cIndex][eIndex][frame] = Channel.CheckError(c,r,gLength);

                    //情報ビットの非誤りビットのインデックス
                    List<Integer> noErrorBitIndex = new ArrayList<>();
                    for(int i = 0;i < g.length;i++){
                        if(c[i] == r[i])noErrorBitIndex.add(i);
                    }

                    //復号時間計測開始時間
                    long startDecode = System.nanoTime();

                    //対数領域sum-product復号,確率領域sum-product復号法,Min-Sum復号法
                    LogDecoder.DecodeResult result = LogDecoder.decode(encodedH,r,eValues[eIndex],maxL);
//                    ProbDecoder.DecodingResult result = ProbDecoder.decode(encodedH,r,eValues[eIndex],maxL);
//                    MinSumDecoder.DecodeResult result = MinSumDecoder.decode(encodedH,r,eValues[eIndex],maxL);

                    //復号時間計測終了時間
                    long endDecode = System.nanoTime();
                    long diffDecodeTime = endDecode - startDecode;

                    //復号後と反復回数とシンドロームの取得
                    int[] decodedC = result.decodedCode();
                    int iterations = result.iterationNum();
                    int syndrome = result.syndrome();

                    //フレームの正誤判定
                    boolean isFrameCorrect = Arrays.equals(c,decodedC); //trueなら成功.falseなら失敗

                    //集計処理
                    synchronized(lock) {
                        //復号時間の加算
                        sumDecodeTime[0] += diffDecodeTime;

                        //反復回数の保存
                        iterationDistribution[cIndex][eIndex][iterations-1] ++;

                        //正誤毎の反復回数,フレーム数の加算
                        if(isFrameCorrect){
                            trueIterations[0] += iterations;
                            trueIterations[1] ++;
                        }else{
                            falseIterations[0] += iterations;
                            falseIterations[1] ++;
                        }

                        //情報ビットの残留誤りビットと誤訂正ビットの加算
                        for(int z = 0;z < g.length;z++){
                            if(c[z] != decodedC[z]){
                                residualsErrorInfoBits[cIndex][eIndex] ++;
                                currentInfoFrameErrorBits++;
                            }
                        }
                        for(int x : noErrorBitIndex){
                            if(c[x] != decodedC[x])errorCorrectionBits[cIndex][eIndex] ++;
                        }

                        //検査ビットの正誤判定
                        for(int z = g.length;z < n;z++){
                            if(c[z] != decodedC[z]){
                                errorCheckBitsCounter[0] ++;
                                currentCheckFrameErrorBits++;
                            }
                        }
                        if(currentCheckFrameErrorBits != 0) errorCheckFrameCounter[0]++;

                        //情報ビットの正誤
                        if(currentInfoFrameErrorBits != 0) errorInfoFrameCounter[0]++;

                        //実際の誤り率の加算
                        sumChannelBitError[cIndex][eIndex] += actualChannelBitErrorRate[cIndex][eIndex][frame];

                        //訂正したと勘違いした数
                        if(syndrome == 0 && !isFrameCorrect) undetectedErrors[cIndex][eIndex]++;
                    }
                }); //並列処理終了

                //復号時間の保存
                decodeTimes[column][errorRate] = sumDecodeTime[0] / 60_000_000_000.0;
                executionTimes[column][1] += decodeTimes[column][errorRate];

                //正誤毎の反復回数の平均
                averageTrueIterations[column][errorRate] =
                        (trueIterations[1] == 0) ? Double.NaN : (double) trueIterations[0] / trueIterations[1];
                averageFalseIterations[column][errorRate] =
                        (falseIterations[1] == 0) ? Double.NaN : (double) falseIterations[0] / falseIterations[1];

                //FER
                frameErrorRate[column][errorRate] = (double)falseIterations[1]/numFrames;

                //IBER
                infoBitErrorRate[column][errorRate] = (double)residualsErrorInfoBits[column][errorRate] / (g.length * numFrames);

                //IFER
                informationFrameErrorRate[column][errorRate] = (double)errorInfoFrameCounter[0]/numFrames;

                //CFER
                checkFrameErrorRate[column][errorRate] = (double)errorCheckFrameCounter[0]/numFrames;

                //平均検査ビット誤り数
                avgErrorCheckBits[column][errorRate] = (double)errorCheckBitsCounter[0]/errorCheckFrameCounter[0];

                //実際の誤り率の平均
                aveChannelBitErrorRate[column][errorRate] = sumChannelBitError[column][errorRate]/numFrames;
            }
            //実行全体の計測終了
            long endTotal = System.currentTimeMillis();
            executionTimes[column][0] = (endTotal - startTotal) / 60_000.0;
            executionTimes[column][2] = (executionTimes[column][1] / executionTimes[column][0]) * 100;
        }

        //実際の通信路誤り率の分散を求める
        for(int i = 0;i < numCM;i++){
            for(int j = 0;j < eValues.length;j++){
                for(int k = 0;k < numFrames;k++){
                    varianceChannelBitError[i][j] += Math.pow((actualChannelBitErrorRate[i][j][k] - aveChannelBitErrorRate[i][j]),2);
                }
                varianceChannelBitError[i][j] /= numFrames;
            }
        }

        //ファイルへの書き出し
        try (PrintWriter pw = new PrintWriter(fileNames, Charset.forName("Windows-31j"))){

            for(int i = 0;i < numCM;i++)pw.printf("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,,,,,,","符号長","行重み","列重み","符号化率","最大反復回数","フレーム数","全体時間(m)","合計復号時間(m)","復号割合(%)","行列番号");
            pw.printf("\n");

            for(int i = 0;i < numCM;i++) pw.printf("%d,%d,%d,%.6f,%d,%d,%.2f,%.2f,%.2f,%d,,,,,,",
                    n,wr,wc,(1-(double)wc/wr),maxL,numFrames,executionTimes[i][0],executionTimes[i][1],executionTimes[i][2],i);
            pw.printf("\n\n");

            for(int i = 0;i < numCM;i++)pw.printf("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,,","通信路誤り率","実際の通信路誤り率の平均","実際の通信路誤り率の分散",
                    "FER","IFER","CFER","平均誤り検査ビット数","s=0だが誤訂正","IBER","成功時の平均繰り返し回数","失敗時の平均繰り返し回数","平均誤訂正ビット率","誤訂正ビット/残留ビット","各誤り率の復号時間(m)");
            pw.printf("\n");

            for(int i = 0;i < eValues.length;i++){
                for(int j = 0;j < numCM;j++){

                    double misRate = (residualsErrorInfoBits[j][i] == 0) ? 0.0 : ((double)errorCorrectionBits[j][i] / residualsErrorInfoBits[j][i]);

                    pw.printf("%.3f,%.6e,%.10e,%.6e,%.6e,%6e,%10e,%d,%.6e,%.4f,%.4f,%.6e,(%d/%d),%.2f,,",
                            eValues[i],
                            aveChannelBitErrorRate[j][i],
                            varianceChannelBitError[j][i],
                            frameErrorRate[j][i],
                            informationFrameErrorRate[j][i],
                            checkFrameErrorRate[j][i],
                            avgErrorCheckBits[j][i],
                            undetectedErrors[j][i],
                            infoBitErrorRate[j][i],
                            averageTrueIterations[j][i],
                            averageFalseIterations[j][i],
                            misRate,
                            errorCorrectionBits[j][i],
                            residualsErrorInfoBits[j][i],
                            decodeTimes[j][i]);
                }
                pw.printf("\n");
            }

            pw.printf("\n以下は反復回数の度数分布\n");
            for(int i = 0;i < numCM;i++){
                pw.printf("回数\\通信路誤り率,");
                for (double eValue : eValues) {
                    pw.printf("%.2f,", eValue);
                }
                pw.printf(",,");
            }
            pw.printf("\n");

            for(int i = 0;i < maxL;i++){
                for(int j = 0;j < numCM;j++){
                    pw.printf("%d,",i);
                    for(int k = 0;k < eValues.length;k++){
                        pw.printf("%d,",iterationDistribution[j][k][i]);
                    }
                    pw.printf(",,");
                }
                pw.printf("\n");
            }
    } catch (IOException e) {
            e.printStackTrace();
        }
    }
}