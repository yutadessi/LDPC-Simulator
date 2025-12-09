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
        String fileNAMEME = "Try";
        //------------------------------
        String fileNames = fileNAMEME + "-result.csv"; //結果保存ファイル名

        //符号パラメータ
        int n = 1024; //符号長
        int wr = 8; //行重み(n % wr = 0)
        int[] wc = {4,5}; //列重み
        int maxL = 50; //最大反復回数
        int numFrames = 100; //フレーム数

        //通信路誤り率eの集合
        double[] eValues = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
//            double[] e = {0.05};

        //タイム計測用配列
        double[][] executionTimes = new double[wc.length][3]; //トータルの実行時間
        double[][] decodeTimes = new double[wc.length][eValues.length]; //各誤り率の復号時間

        //出力用保存配列-各誤り率のデータ
        double[][][] actualChannelBitErrorRate = new double[wc.length][eValues.length][numFrames]; //各フレームの実際の通信路誤り率
        double[][] sumChannelBitError = new double[wc.length][eValues.length]; //各誤り率における実際の通信路誤り率の合計
        double[][] aveChannelBitErrorRate = new double[wc.length][eValues.length]; //各誤り率における実際の通信路誤り率の平均
        double[][] varianceChannelBitError = new double[wc.length][eValues.length]; //各誤り率における実際の通信路誤り率の分散
        double[][] frameErrorRate = new double[wc.length][eValues.length]; //FER(符号長の失敗確率)
        double[][] infomationFrameErrorRate = new double[wc.length][eValues.length]; //IFER(情報ビットのみの失敗率率）
        double[][] infoBitErrorRate = new double[wc.length][eValues.length]; //IBER(Info Bit Error Rate)
        double[][] averageTrueIterations = new double[wc.length][eValues.length]; //訂正成功時の平均繰り返し回数
        double[][] averageFalseIterations = new double[wc.length][eValues.length]; //訂正失敗時の平均繰り返し回数
        int[][] residualsErrorBits = new int[wc.length][eValues.length]; //情報ビットの残留誤りビット数
        int[][] errorCorrectionBits = new int[wc.length][eValues.length]; //情報ビットの誤訂正ビット数
        double[][][] iterationDistribution = new double[wc.length][eValues.length][maxL]; //反復回数の度数分布
        int[][] undetectedErrors = new int[wc.length][eValues.length]; //シンドロームは0だが,誤訂正している数

        //各列重みでのシミュレーション実行
        for(int column = 0;column < wc.length;column++){

            //実行時間全体の計測開始
            long startTotal = System.currentTimeMillis();

            //検査行列Hと生成行列Gの作成
            int [][] h = GenerateMatrix.gallagerCheckMatrix(n,wr,wc[column]);
            int [][] g = GenerateMatrix.generatorMatrix(h);

            //検査行列を保存
            String filePath = fileNAMEME + wc[column] + "-HMatrix.txt"; //検査行列保存ファイル名
            CheckMatrixIO.saveCheckMatrix(h,filePath);

            //HとGを組織符号化
            List<Integer> columnIndicatesToSwap = new ArrayList<>();  //組織符号化用インデックス
            int[][] encodedG = TissueEncoder.EncodeG(g,columnIndicatesToSwap);
            int[][] encodedH = TissueEncoder.EncodeH(h,columnIndicatesToSwap);

            //各通信路誤り率でのシミュレーション
            for(int errorRate = 0; errorRate < eValues.length; errorRate++){

                //復号時間の合計用変数を初期化
                long sumDecodeTime = 0;

                //正誤毎の反復回数,フレーム数の合計([0]は反復回数,[1]はフレーム数)
                int[] trueIterations = new int[2];
                int[] falseIterations = new int[2];

                int errorInfoBitsCounter = 0;

                for(int frame = 0;frame < numFrames;frame++){
                     //メッセージと送信語、受信語の作成
                    int[] c = GenerateC.geneC(encodedG);
                    int[] r = Channel.GenerateR(c,eValues[errorRate]);

                    //フレームごとの情報ビットの正誤
                    int currentInfoFrameErrorBits = 0;

                    //実際の通信路での誤り率の取得
                    actualChannelBitErrorRate[column][errorRate][frame] = Channel.CheckError(c,r);

                    //情報ビットの非誤りビットのインデックス
                    List<Integer> noErrorBitIndex = new ArrayList<>();
                    for(int i = 0;i < g.length;i++){
                        if(c[i] == r[i])noErrorBitIndex.add(i);
                    }

                    //復号時間計測開始時間
                    long startDecode = System.nanoTime();

                    //対数領域sum-product復号,確率領域sum-product復号法,Min-Sum復号法
                    LogDecoder.DecodeResult result = LogDecoder.decode(encodedH,r,eValues[errorRate],maxL);
//                    ProbDecoder.DecodingResult result = ProbDecoder.decode(encodedH,r,eValues[errorRate],maxL);
//                    MinSumDecoder.DecodeResult result = MinSumDecoder.decode(encodedH,r,eValues[errorRate],maxL);

                    //復号時間計測終了時間と復号時間の加算
                    long endDecode = System.nanoTime();
                    sumDecodeTime += (endDecode - startDecode);

                    //復号後と反復回数とシンドロームの取得
                    int[] decodedC = result.decodedCode();
                    int iterations = result.iterationNum();
                    int syndrome = result.syndrome();

                    //フレームの正誤判定
                    boolean isFrameCorrect = Arrays.equals(c,decodedC); //trueなら成功.falseなら失敗

                    //反復回数の保存
                    iterationDistribution[column][errorRate][iterations-1] ++;

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
                            residualsErrorBits[column][errorRate] ++;
                            currentInfoFrameErrorBits++;
                        }
                    }
                    for(int x : noErrorBitIndex){
                        if(c[x] != decodedC[x])errorCorrectionBits[column][errorRate] ++;
                    }

                    //情報ビットの正誤
                    if(currentInfoFrameErrorBits != 0)errorInfoBitsCounter++;

                    //実際の誤り率の加算
                    sumChannelBitError[column][errorRate] += actualChannelBitErrorRate[column][errorRate][frame];

                    //訂正したと勘違いした数
                    if(syndrome == 0 && !isFrameCorrect) undetectedErrors[column][errorRate]++;
                }

                //復号時間の保存
                decodeTimes[column][errorRate] = sumDecodeTime / 60_000_000_000.0;
                executionTimes[column][1] += decodeTimes[column][errorRate];

                //正誤毎の反復回数の平均
                averageTrueIterations[column][errorRate] = (double)trueIterations[0]/ trueIterations[1];
                averageFalseIterations[column][errorRate] = (double)falseIterations[0]/ falseIterations[1];

                //FER
                frameErrorRate[column][errorRate] = (double)falseIterations[1]/numFrames;

                //IBER
                infoBitErrorRate[column][errorRate] = (double)residualsErrorBits[column][errorRate] / (g.length * numFrames);

                //IFER
                infomationFrameErrorRate[column][errorRate] = (double)errorInfoBitsCounter/numFrames;

                //実際の誤り率の平均
                aveChannelBitErrorRate[column][errorRate] = sumChannelBitError[column][errorRate]/numFrames;
            }
            //実行全体の計測終了(書き出しを除く)
            long endTotal = System.currentTimeMillis();
            executionTimes[column][0] = (endTotal - startTotal) / 60_000.0;
            executionTimes[column][2] = (executionTimes[column][1] / executionTimes[column][0]) * 100;
        }

        //実際の通信路誤り率の分散を求める
        for(int i = 0;i < wc.length;i++){
            for(int j = 0;j < eValues.length;j++){
                for(int k = 0;k < numFrames;k++){
                    varianceChannelBitError[i][j] += Math.pow((actualChannelBitErrorRate[i][j][k] - aveChannelBitErrorRate[i][j]),2);
                }
                varianceChannelBitError[i][j] /= numFrames;
            }
        }

        long startWrite = System.nanoTime();

        try (PrintWriter pw = new PrintWriter(fileNames, Charset.forName("Windows-31j"))){
            for(int i = 0;i < wc.length;i++)pw.printf("%s,%s,%s,%s,%s,%s,%s,%s,%s,,,,|,","符号長","行重み","列重み","符号化率","最大反復回数","フレーム数","全体時間(m)","合計復号時間(m)","復号割合(%)");
            pw.printf("\n");
            for(int i = 0;i < wc.length;i++) pw.printf("%s,%s,%s,%s,%s,%s,%.2f,%.2f,%.2f,,,,|,",n,wr,wc[i],(1-(double)wc[i]/wr),maxL,numFrames,executionTimes[i][0],executionTimes[i][1],executionTimes[i][2]);
            pw.printf("\n\n");
            for(int i = 0;i < wc.length;i++)pw.printf("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,|,","通信路誤り率","実際の通信路誤り率の平均","実際の通信路誤り率の分散",
                    "FER","IFER","s=0だが誤訂正","IBER","成功時の平均繰り返し回数","失敗時の平均繰り返し回数","平均誤訂正ビット率","誤訂正ビット/残留ビット","各誤り率の復号時間(m)");
            pw.printf("\n");
            for(int i = 0;i < eValues.length;i++){
                for(int j = 0;j < wc.length;j++){
                    pw.printf("%.2f,%s,%.10f,%.4f,%.4f,%s,%s,%s,%s,%.6f,(%s/%s),%.2f,|,", eValues[i], aveChannelBitErrorRate[j][i],varianceChannelBitError[j][i],
                            frameErrorRate[j][i],infomationFrameErrorRate[j][i],undetectedErrors[j][i], infoBitErrorRate[j][i],averageTrueIterations[j][i],averageFalseIterations[j][i],
                            ((double)errorCorrectionBits[j][i]/residualsErrorBits[j][i]),errorCorrectionBits[j][i],residualsErrorBits[j][i],decodeTimes[j][i]);
                }
                pw.printf("\n");
            }
            pw.printf("\n以下は反復回数の度数分布\n");
            for(int i = 0;i < wc.length;i++){
                pw.printf("回数\\通信路誤り率,");
                for(int j = 0;j < eValues.length;j++){
                    pw.printf("%s,",eValues[j]);
                }
                pw.printf(",");
            }
            pw.printf("\n");
            for(int i = 0;i < maxL;i++){
                for(int j = 0;j < wc.length;j++){
                    pw.printf("%s,",i);
                    for(int k = 0;k < eValues.length;k++){
                        pw.printf("%s,",iterationDistribution[j][k][i]);
                    }
                    pw.printf(",");
                }
                pw.printf("\n");
            }
            long endWrite = System.nanoTime();
            pw.printf("書き出し時間（秒）: %s",(endWrite - startWrite) / 1_000_000_000.0);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
