package LDPC;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.io.PrintWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.stream.IntStream;

public class LDPC_LoadIntStream {
    public static void main(String[] args) {

        //ファイル名、毎回変える！！--------
        String fileNAMEME = "No.3";
        //------------------------------

        String fileNames = fileNAMEME + "-LoadHResult.txt";
        String filePath = fileNAMEME + "-HMatrix.txt";

        //符号パラメーター
        int maxL = 20; //最大反復回数
        int numFrames = 10_000; //フレーム数

        //通信路誤り率eの設定
        double[] eValues = {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
//        double[] eValues = {0.03};

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

        //タイム計測用配列
        double[] executionTimes = new double[3]; //トータルの実行時間
        double[] decodeTimes = new double[eValues.length]; //各誤り率の復号時間

        //出力用保存配列-各誤り率のデータ
        double[][] actualChannelBitErrorRate = new double[eValues.length][numFrames]; //各フレームの実際の通信路誤り率
        double[] sumChannelBitError = new double[eValues.length]; //各誤り率における実際の通信路誤り率の合計
        double[] aveChannelBitErrorRate = new double[eValues.length]; //各誤り率における実際の通信路誤り率の平均
        double[] varianceChannelBitError = new double[eValues.length]; //各誤り率における実際の通信路誤り率の分散
        double[] frameErrorRate = new double[eValues.length]; //FER(符号長の失敗確率)
        double[] informationFrameErrorRate = new double[eValues.length]; //IFER(情報ビットのみの失敗率率）
        double[] infoBitErrorRate = new double[eValues.length]; //IBER(Info Bit Error Rate)
        double[] averageTrueIterations = new double[eValues.length]; //訂正成功時の平均繰り返し回数
        double[] averageFalseIterations = new double[eValues.length]; //訂正失敗時の平均繰り返し回数
        int[] residualsErrorBits = new int[eValues.length]; //情報ビットの残留誤りビット数
        int[] errorCorrectionBits = new int[eValues.length]; //情報ビットの誤訂正ビット数
        double[][] iterationDistribution = new double[eValues.length][maxL]; //反復回数の度数分布
        int[] undetectedErrors = new int[eValues.length]; //シンドロームは0だが,誤訂正している数


        //実行時間全体の計測開始
        long startTotal = System.currentTimeMillis();

        //検査行列hの復元と生成行列gの作成
        int [][] h = CheckMatrixIO.loadCheckMatrix(filePath);
        int [][] g = GenerateMatrix.generatorMatrix(h);

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
                actualChannelBitErrorRate[eIndex][frame] = Channel.CheckError(c,r);

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
                    iterationDistribution[eIndex][iterations-1] ++;

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
                            residualsErrorBits[eIndex] ++;
                            currentInfoFrameErrorBits++;
                        }
                    }
                    for(int x : noErrorBitIndex){
                        if(c[x] != decodedC[x])errorCorrectionBits[eIndex] ++;
                    }

                    //情報ビットの正誤
                    if(currentInfoFrameErrorBits != 0) errorInfoBitsCounter[0]++;

                    //実際の誤り率の加算
                    sumChannelBitError[eIndex] += actualChannelBitErrorRate[eIndex][frame];

                    //訂正したと勘違いした数
                    if(syndrome == 0 && !isFrameCorrect) undetectedErrors[eIndex]++;
                }
            }); //並列処理終了

            //復号時間の保存
            decodeTimes[errorRate] = sumDecodeTime[0] / 60_000_000_000.0;
            executionTimes[1] += decodeTimes[errorRate];

            //正誤毎の反復回数の平均
            averageTrueIterations[errorRate] = (double)trueIterations[0]/ trueIterations[1];
            averageFalseIterations[errorRate] = (double)falseIterations[0]/ falseIterations[1];

            //FER
            frameErrorRate[errorRate] = (double)falseIterations[1]/numFrames;

            //IBER
            infoBitErrorRate[errorRate] = (double)residualsErrorBits[errorRate] / (g.length * numFrames);

            //IFER
            informationFrameErrorRate[errorRate] = (double)errorInfoBitsCounter[0]/numFrames;

            //実際の誤り率の平均
            aveChannelBitErrorRate[errorRate] = sumChannelBitError[errorRate]/numFrames;
        }
        //実行全体の計測終了
        long endTotal = System.currentTimeMillis();
        executionTimes[0] = (endTotal - startTotal) / 60_000.0;
        executionTimes[2] = (executionTimes[1] / executionTimes[0]) * 100;


        //実際の通信路誤り率の分散を求める
        for(int j = 0;j < eValues.length;j++){
            for(int k = 0;k < numFrames;k++){
                varianceChannelBitError[j] += Math.pow((actualChannelBitErrorRate[j][k] - aveChannelBitErrorRate[j]),2);
            }
            varianceChannelBitError[j] /= numFrames;
        }


        //ファイルへの書き出し
        try (PrintWriter pw = new PrintWriter(fileNames, Charset.forName("Windows-31j"))){
            pw.printf("%s,%s,%s,%s,%s,%s,%s,%s,%s\n","符号長","行重み","列重み","符号化率","最大反復回数","フレーム数","全体時間(m)","合計復号時間(m)","復号割合(%)");
            pw.printf("%s,%s,%s,%s,%s,%s,%.2f,%.2f,%.2f\n\n",n,wr,wc,(1-(double)wc/wr),maxL,numFrames,executionTimes[0],executionTimes[1],executionTimes[2]);
            pw.printf("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","通信路誤り率","実際の通信路誤り率の平均","実際の通信路誤り率の分散",
                    "FER","IFER","s=0だが誤訂正","IBER","成功時の平均繰り返し回数","失敗時の平均繰り返し回数","平均誤訂正ビット率","誤訂正ビット/残留ビット","各誤り率の復号時間(m)");
            for(int i = 0;i < eValues.length;i++){
                    pw.printf("%.2f,%s,%.10f,%.4f,%.4f,%s,%s,%s,%s,%.6f,(%s/%s),%.2f,,", eValues[i], aveChannelBitErrorRate[i],varianceChannelBitError[i],
                            frameErrorRate[i],informationFrameErrorRate[i],undetectedErrors[i], infoBitErrorRate[i],averageTrueIterations[i],averageFalseIterations[i],
                            ((double)errorCorrectionBits[i]/residualsErrorBits[i]),errorCorrectionBits[i],residualsErrorBits[i],decodeTimes[i]);
                pw.printf("\n");
            }
            pw.printf("\n以下は反復回数の度数分布\n");
            pw.printf("回数\\通信路誤り率,");
            for (double eValue : eValues) pw.printf("%s,", eValue);
            pw.printf("\n");
            for(int i = 0;i < maxL;i++){
                pw.printf("%s,",i);
                for(int k = 0;k < eValues.length;k++)pw.printf("%s,",iterationDistribution[k][i]);
                pw.printf("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}