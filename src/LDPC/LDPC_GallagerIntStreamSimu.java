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
        String fileNAMEME = "8-4(10_000)Gallager-Log";
        //------------------------------
        String fileNames = fileNAMEME + "-result.csv";

        //符号パラメータ
        int n = 1024; //符号長
        int wr = 8; //行重み(n % wr = 0)
        int wc = 4; //列重み
        int maxL = 50; //最大反復回数
        int numFrames = 10; //フレーム数

        //検査行列作成数
        int numCM = 50;

        //通信路誤り率eの集合
        double[] eValues = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
//            double[] e = {0.05};

        //タイム計測用配列
        double[][] executionTimes = new double[numCM][3]; //トータルの実行時間
        double[][] decodeTimes = new double[numCM][eValues.length]; //各誤り率の復号時間

        //出力用保存配列-各誤り率のデータ
        double[][][] actualChannelBitErrorRate = new double[numCM][eValues.length][numFrames]; //各フレームの実際の通信路誤り率
        double[][] sumChannelBitError = new double[numCM][eValues.length]; //各誤り率における実際の通信路誤り率の合計
        double[][] aveChannelBitErrorRate = new double[numCM][eValues.length]; //各誤り率における実際の通信路誤り率の平均
        double[][] varianceChannelBitError = new double[numCM][eValues.length]; //各誤り率における実際の通信路誤り率の分散
        double[][] frameErrorRate = new double[numCM][eValues.length]; //FER(符号長の失敗確率)
        double[][] informationFrameErrorRate = new double[numCM][eValues.length]; //IFER(情報ビットのみの失敗率率）
        double[][] infoBitErrorRate = new double[numCM][eValues.length]; //IBER(Info Bit Error Rate)
        double[][] averageTrueIterations = new double[numCM][eValues.length]; //訂正成功時の平均繰り返し回数
        double[][] averageFalseIterations = new double[numCM][eValues.length]; //訂正失敗時の平均繰り返し回数
        int[][] residualsErrorBits = new int[numCM][eValues.length]; //情報ビットの残留誤りビット数
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

                final int[] errorInfoBitsCounter = {0};

                Object lock = new Object();

                //並列処理開始
                IntStream.range(0, numFrames).parallel().forEach(frame -> {

                    //メッセージと送信語、受信語の作成
                    int[] c = GenerateC.geneC(encodedG);
                    int[] r = Channel.GenerateR(c,eValues[eIndex]);

                    //フレームごとの情報ビットの正誤
                    int currentInfoFrameErrorBits = 0;

                    //実際の通信路での誤り率の取得
                    actualChannelBitErrorRate[cIndex][eIndex][frame] = Channel.CheckError(c,r);

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

                        //残留誤りビットと誤訂正ビットの加算
                        for(int z = 0;z < g.length;z++){
                            if(c[z] != decodedC[z]){
                                residualsErrorBits[cIndex][eIndex] ++;
                                currentInfoFrameErrorBits++;
                            }
                        }
                        for(int x : noErrorBitIndex){
                            if(c[x] != decodedC[x])errorCorrectionBits[cIndex][eIndex] ++;
                        }

                        //情報ビットの正誤
                        if(currentInfoFrameErrorBits != 0) errorInfoBitsCounter[0]++;

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
                infoBitErrorRate[column][errorRate] = (double)residualsErrorBits[column][errorRate] / (g.length * numFrames);

                //IFER
                informationFrameErrorRate[column][errorRate] = (double)errorInfoBitsCounter[0]/numFrames;

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

        // ===== 最終スコア計算（小さいほど良い） =====
        final double eps = 1e-15;

        // eValues方向に集約した特徴量（cmごと）
        double[] avgLogFER_CM = new double[numCM];
        double[] avgLogIBER_CM = new double[numCM];
        double[] avgIter_CM = new double[numCM];
        double[] avgTime_CM = new double[numCM];

        for (int cm = 0; cm < numCM; cm++) {
            double sLogFer = 0.0, sLogIber = 0.0, sIter = 0.0, sTime = 0.0;

            for (int ei = 0; ei < eValues.length; ei++) {
                sLogFer += Math.log10(frameErrorRate[cm][ei] + eps);
                sLogIber += Math.log10(infoBitErrorRate[cm][ei] + eps);
                // 全フレームの平均反復回数（成功/失敗に依存せず必ず定義できる）
                int iterSum = 0;
                for (int it = 0; it < maxL; it++) {
                    iterSum += (it + 1) * iterationDistribution[cm][ei][it]; // it=0 が 1回
                }
                double avgIterAll = (double) iterSum / numFrames;
                sIter += avgIterAll;

                sTime += decodeTimes[cm][ei];
            }

            int m = eValues.length;
            avgLogFER_CM[cm] = sLogFer / m;
            avgLogIBER_CM[cm] = sLogIber / m;
            avgIter_CM[cm] = sIter / m;
            avgTime_CM[cm] = sTime / m;
        }

        // min-max 正規化（0..1）
        double minFer = Double.POSITIVE_INFINITY, maxFer = Double.NEGATIVE_INFINITY;
        double minIber = Double.POSITIVE_INFINITY, maxIber = Double.NEGATIVE_INFINITY;
        double minIter = Double.POSITIVE_INFINITY, maxIter = Double.NEGATIVE_INFINITY;
        double minTime = Double.POSITIVE_INFINITY, maxTime = Double.NEGATIVE_INFINITY;

        for (int cm = 0; cm < numCM; cm++) {
            if (avgLogFER_CM[cm] < minFer) minFer = avgLogFER_CM[cm];
            if (avgLogFER_CM[cm] > maxFer) maxFer = avgLogFER_CM[cm];

            if (avgLogIBER_CM[cm] < minIber) minIber = avgLogIBER_CM[cm];
            if (avgLogIBER_CM[cm] > maxIber) maxIber = avgLogIBER_CM[cm];

            if (avgIter_CM[cm] < minIter) minIter = avgIter_CM[cm];
            if (avgIter_CM[cm] > maxIter) maxIter = avgIter_CM[cm];

            if (avgTime_CM[cm] < minTime) minTime = avgTime_CM[cm];
            if (avgTime_CM[cm] > maxTime) maxTime = avgTime_CM[cm];
        }

        double rangeFer = maxFer - minFer;
        double rangeIber = maxIber - minIber;
        double rangeIter = maxIter - minIter;
        double rangeTime = maxTime - minTime;

        double[] score_CM = new double[numCM];

        for (int cm = 0; cm < numCM; cm++) {
            double nFer  = (rangeFer  == 0.0) ? 0.0 : (avgLogFER_CM[cm]  - minFer)  / rangeFer;
            double nIber = (rangeIber == 0.0) ? 0.0 : (avgLogIBER_CM[cm] - minIber) / rangeIber;
            double nIter = (rangeIter == 0.0) ? 0.0 : (avgIter_CM[cm]    - minIter) / rangeIter;
            double nTime = (rangeTime == 0.0) ? 0.0 : (avgTime_CM[cm]    - minTime) / rangeTime;

            // 重み（必要なら調整）
            score_CM[cm] = 0.55 * nFer + 0.25 * nIber + 0.10 * nIter + 0.10 * nTime;
        }

        // Score の統計（平均/最小/最大/分散）と最良行列番号
        double sumScore = 0.0;
        double minScore = Double.POSITIVE_INFINITY;
        double maxScore = Double.NEGATIVE_INFINITY;
        int bestCM = -1;

        for (int cm = 0; cm < numCM; cm++) {
            double s = score_CM[cm];
            sumScore += s;
            if (s < minScore) { minScore = s; bestCM = cm; }
            if (s > maxScore) maxScore = s;
        }
        double meanScore = sumScore / numCM;

        double varScore = 0.0;
        for (int cm = 0; cm < numCM; cm++) {
            double d = score_CM[cm] - meanScore;
            varScore += d * d;
        }
        varScore /= numCM; // 母分散



        //ファイルへの書き出し
        try (PrintWriter pw = new PrintWriter(fileNames, Charset.forName("Windows-31j"))){

            for(int i = 0;i < numCM;i++)pw.printf("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,,,,","符号長","行重み","列重み","符号化率","最大反復回数","フレーム数","全体時間(m)","合計復号時間(m)","復号割合(%)","行列番号");
            pw.printf("\n");

            for(int i = 0;i < numCM;i++) pw.printf("%d,%d,%d,%.6f,%d,%d,%.2f,%.2f,%.2f,%d,,,,",
                    n,wr,wc,(1-(double)wc/wr),maxL,numFrames,executionTimes[i][0],executionTimes[i][1],executionTimes[i][2],i);
            pw.printf("\n\n");

            for(int i = 0;i < numCM;i++)pw.printf("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,,","通信路誤り率","実際の通信路誤り率の平均","実際の通信路誤り率の分散",
                    "FER","IFER","s=0だが誤訂正","IBER","成功時の平均繰り返し回数","失敗時の平均繰り返し回数","平均誤訂正ビット率","誤訂正ビット/残留ビット","各誤り率の復号時間(m)");
            pw.printf("\n");

            for(int i = 0;i < eValues.length;i++){
                for(int j = 0;j < numCM;j++){

                    double misRate = (residualsErrorBits[j][i] == 0) ? 0.0 : ((double)errorCorrectionBits[j][i] / residualsErrorBits[j][i]);

                    pw.printf("%.2f,%.6e,%.10e,%.6e,%.6e,%d,%.6e,%.4f,%.4f,%.6e,(%d/%d),%.2f,,",
                            eValues[i],
                            aveChannelBitErrorRate[j][i],
                            varianceChannelBitError[j][i],
                            frameErrorRate[j][i],
                            informationFrameErrorRate[j][i],
                            undetectedErrors[j][i],
                            infoBitErrorRate[j][i],
                            averageTrueIterations[j][i],
                            averageFalseIterations[j][i],
                            misRate,
                            errorCorrectionBits[j][i],
                            residualsErrorBits[j][i],
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
                        pw.printf("%d,",(int)iterationDistribution[j][k][i]); // int化してなくても動くようにキャスト
                    }
                    pw.printf(",,");
                }
                pw.printf("\n");
            }

            // ===== 最終スコア（追記） =====
            pw.printf("\n最終スコア（小さいほど良い）\n");
            pw.printf("行列番号,Score,AvgLogFER,AvgLogIBER,AvgIter,AvgDecodeTime(m)\n");
            for (int cm = 0; cm < numCM; cm++) {
                pw.printf("%d,%.6f,%.6f,%.6f,%.4f,%.6f\n",
                        cm, score_CM[cm], avgLogFER_CM[cm], avgLogIBER_CM[cm], avgIter_CM[cm], avgTime_CM[cm]);
            }
            pw.printf("Score統計,mean=%.6f,min=%.6f,max=%.6f,var=%.6f,bestCM=%d\n",
                    meanScore, minScore, maxScore, varScore, bestCM);
    } catch (IOException e) {
            e.printStackTrace();
        }
    }
    static double mean(double[] a) {
        double s = 0;
        for (double v : a) s += v;
        return s / a.length;
    }

    static double variance(double[] a) {
        double m = mean(a);
        double s = 0;
        for (double v : a) {
            double d = v - m;
            s += d * d;
        }
        return s / a.length; // 母分散
    }

    static double min(double[] a) {
        double m = Double.POSITIVE_INFINITY;
        for (double v : a) if (v < m) m = v;
        return m;
    }

    static double max(double[] a) {
        double m = Double.NEGATIVE_INFINITY;
        for (double v : a) if (v > m) m = v;
        return m;
    }

}