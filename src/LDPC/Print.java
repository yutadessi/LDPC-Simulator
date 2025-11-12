package LDPC;

public class Print {
    public static void Matrix(int[][] Mat){
        System.out.println("Matrix [" + Mat.length + ", " + Mat[0].length + "]");
        for(int i = 0;i < Mat.length;i++){
            for(int j = 0;j < Mat[0].length;j++) {
                System.out.print(Mat[i][j] + " ");
            }
            System.out.println(" ");
        }
        System.out.println(" ");
    }

    public static void Array(int[] Array){
        System.out.println("Array [" + Array.length + "]");
        for(int i = 0;i < Array.length;i++){
            System.out.print(Array[i] + " ");
        }
        System.out.println(" \n");
    }
}
