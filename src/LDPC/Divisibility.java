package LDPC;

public class Divisibility {
    public static int checker(int n,int wr) {
        if(n % wr == 0)return 0;
        System.out.println("符号長が行重みで割り切れません");
        return 1;
    }
}
