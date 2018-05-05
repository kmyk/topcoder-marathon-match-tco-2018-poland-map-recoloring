import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.security.*;
import java.util.*;
import javax.imageio.*;

public class MapRecoloringVis {
    static int minS = 20, maxS = 200;
    static int minC = 2, maxC = 5;
    static int[] dr = {1, 0,-1, 0};
    static int[] dc = {0, 1, 0,-1};

    int H, W;               // size of the map
    int R;                  // number of regions
    int[][] regions;        // which region this pixel of the map belongs to
    int[][] oldColors;      // old (starting) colors of the regions of the map (one color per pixel)
    int[] newColors;        // colors of the regions assigned by the solution (one color per region)

    int[] regPrimaryColor; // primary colors of the regions (for visualization purposes)
    // -----------------------------------------
    boolean isInside(int r, int c) {
        return (r >= 0 && r < H && c >= 0 && c < W);
    }
    // -----------------------------------------
    int getColor(boolean isNew, int r, int c) {
        return isNew ? newColors[regions[r][c]] : oldColors[r][c];
    }
    // -----------------------------------------
    String generate(String seedStr) {
      try {
        // generate test case
        SecureRandom rnd = SecureRandom.getInstance("SHA1PRNG"); 
        long seed = Long.parseLong(seedStr);
        rnd.setSeed(seed);
        H = rnd.nextInt(maxS - minS + 1) + minS;
        W = rnd.nextInt(maxS - minS + 1) + minS;
        int C = rnd.nextInt(maxC - minC + 1) + minC;
        if (seed == 1) {
            W = H = minS;
            C = 2;
        }
        else if (seed == 2) {
            W = H = maxS;
            C = maxC;
        }

        // pick the number of regions based on the total # of pixels in the map
        // the average size of the region will be between 10 and 50 pixels
        // max # of regions can be W * H / 10 = 4k
        int szr = rnd.nextInt(41) + 10;
        R = W * H / szr;

        // generate the map of regions
        regions = new int[H][W];
        for (int i = 0; i < H; ++i)
            Arrays.fill(regions[i], -1);

        // seed the regions (some regions can have several seeds and end up disconnected)
        for (int i = 0; i < R; ++i) {
            int parts = rnd.nextInt(4) + 1;
            for (int j = 0; j < parts; ++j) {
                int r, c;
                do {
                    r = rnd.nextInt(H - 2) + 1;
                    c = rnd.nextInt(W - 2) + 1;
                } while (regions[r][c] != -1);
                regions[r][c] = i;
            }
        }

        // grow the regions from their seeds
        java.util.List<Integer> candidates = new ArrayList<Integer>();
        for (int r = 0; r < H; ++r)
        for (int c = 0; c < W; ++c) 
            if (regions[r][c] != -1) {
                // unassigned neighbors of this pixel are initial growth candidates
                for (int k = 0; k < 4; ++k) {
                    int nr = r + dr[k], nc = c + dc[k];
                    if (isInside(nr, nc) && regions[nr][nc] == -1 && !candidates.contains(nr * W + nc))
                        candidates.add(nr * W + nc);
                }
            }
        while (!candidates.isEmpty()) {
            // pick a random pixel and assign it to the region of one of its neighbors at random
            int ind = rnd.nextInt(candidates.size());
            int r = candidates.get(ind) / W, c = candidates.get(ind) % W;
            java.util.List<Integer> adjRegions = new ArrayList<Integer>();
            for (int k = 0; k < 4; ++k) {
                int nr = r + dr[k], nc = c + dc[k];
                if (isInside(nr, nc) && regions[nr][nc] != -1)
                    adjRegions.add(k);
            }
            int adjK = adjRegions.get(rnd.nextInt(adjRegions.size()));
            regions[r][c] = regions[r + dr[adjK]][c + dc[adjK]];
            // remove this pixel from the candidates
            candidates.remove(ind);
            // and add its empty neighbors to the candidates list if they are not there yet
            for (int k = 0; k < 4; ++k) {
                int nr = r + dr[k], nc = c + dc[k];
                if (isInside(nr, nc) && regions[nr][nc] == -1 && !candidates.contains(nr * W + nc))
                    candidates.add(nr * W + nc);
            }
        }

        while (true) {
            // color each region into one primary color and several secondary
            regPrimaryColor = new int[R];
            int[][] regSecondaryColors = new int[R][];
            int[] solid = new int[R];
            for (int i = 0; i < R; ++i) {
                regPrimaryColor[i] = rnd.nextInt(C);
                int nSec = rnd.nextInt(C);
                regSecondaryColors[i] = new int[nSec];
                for (int j = 0; j < nSec; ++j) {
                    do { regSecondaryColors[i][j] = rnd.nextInt(C); }
                    while (regSecondaryColors[i][j] == regPrimaryColor[i]);
                }
                solid[i] = (nSec == 0 ? 100 : (rnd.nextInt(50) + 50));
            }

            // assign color to each pixel of each region
            oldColors = new int[H][W];
            for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                int reg = regions[i][j];
                if (rnd.nextInt(100) < solid[reg])
                    oldColors[i][j] = regPrimaryColor[reg];
                else
                    oldColors[i][j] = regSecondaryColors[reg][rnd.nextInt(regSecondaryColors[reg].length)];
            }

            // make sure that the initial map coloring is invalid, i.e. recoloring 0 pixels will not produce a solution
            // check that a pair of adjacent regions with the same-colored pixels on the border exists
            // (if the borders of all regions are ok but some regions are colored in more than one color, regenerate anyways)
            int r1, r2;
            boolean ok = false;
            for (int i = 0; i < H && !ok; ++i)
            for (int j = 0; j < W && !ok; ++j) {
                if (i > 0) {
                    r1 = regions[i][j];
                    r2 = regions[i-1][j];
                    if (r1 != r2 && getColor(false, i, j) == getColor(false, i-1, j)) {
                        ok = true;
                    }
                }
                if (j > 0) {
                    r1 = regions[i][j];
                    r2 = regions[i][j-1];
                    if (r1 != r2 && getColor(false, i, j) == getColor(false, i, j-1)) {
                        ok = true;
                    }
                }
            }
            if (ok)
                break;
        }

        if (vis) {
            // generate palette for R colors; first few are defined and the rest are random
            palette = new ArrayList<Color>();
            int[] firstColors = {0x6633ff, 0xcc33ff, 0xffcc33, 0x33ccff, 0x66ff33};
            for (int i = 0; i < R && i < firstColors.length; ++i)
                palette.add(new Color(firstColors[i]));
            for (int i = firstColors.length; i < R; ++i)
                palette.add(new Color(rnd.nextInt(0x1000000)));
        }

        StringBuffer sb = new StringBuffer();
        sb.append("H = ").append(H).append('\n');
        sb.append("W = ").append(W).append('\n');
        sb.append("R = ").append(R).append('\n');
        sb.append("Regions:\n");
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j)
                sb.append(regions[i][j]).append(" ");
            sb.append('\n');
        }
        sb.append("Old map colors:\n");
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j)
                sb.append(oldColors[i][j]).append(" ");
            sb.append('\n');
        }
        return sb.toString();
      }
      catch (Exception e) {
        addFatalError("An exception occurred while generating test case.");
        e.printStackTrace(); 
        return "";
      }
    }
    // -----------------------------------------
    public double runTest(String seed) {
      try {
        String test = generate(seed);
        if (debug)
            System.out.println(test);

        if (vis) {
            SZX = W * SZ + 1;
            SZY = H * SZ + 1;
            // skip if exists
            if (! (new File("vis/" + fileName + "-in.png").exists())) {
                // draw the starting data
                draw(false);
            }
        }

        if (proc != null) {
            // call the solution
            int[] regionsArg = new int[H * W];
            int[] oldColorsArg = new int[H * W];
            for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                regionsArg[i * W + j] = regions[i][j];
                oldColorsArg[i * W + j] = oldColors[i][j];
            }

            try { 
                newColors = recolor(H, regionsArg, oldColorsArg);
            } catch (Exception e) {
                addFatalError("Failed to get result from recolor.");
                return -1.0;
            }

            // check the return and score it
            if (newColors == null) {
                addFatalError("Your return contained invalid number of elements.");
                return -1.0;
            }
            if (newColors.length != R) {
                addFatalError("Your return contained " + newColors.length + " elements, and it should have contained " + R + ".");
                return -1.0;
            }

            // you never need more than R different colors to color the map
            // so each color can be restricted to 0..R-1
            for (int i = 0; i < R; ++i) {
                if (newColors[i] < 0 || newColors[i] >= R) {
                    addFatalError("Each color in your return must be between 0 and " + (R-1) + ", inclusive.");
                    return -1.0;
                }
            }
        } else {
            // to show off visualization, use the original assignment of primary colors
            newColors = new int[R];
            for (int i = 0; i < R; ++i)
                newColors[i] = regPrimaryColor[i];
        }

        if (vis) {
            // draw the results of program execution (even with invalid score)
            draw(true);
        }

        // check that all pairs of adjacent regions have different colors
        int r1, r2;
        for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            if (i > 0) {
                r1 = regions[i][j];
                r2 = regions[i-1][j];
                if (r1 != r2 && newColors[r1] == newColors[r2]) {
                    addFatalError("Regions "+r1+" and "+r2+ " adjacent in pixels ("+(i-1)+","+j+") and ("+i+","+j+") are of the same color " + newColors[r1] + ".");
                    return -1.0;
                }
            }
            if (j > 0) {
                r1 = regions[i][j];
                r2 = regions[i][j-1];
                if (r1 != r2 && newColors[r1] == newColors[r2]) {
                    addFatalError("Regions "+r1+" and "+r2+ " adjacent in pixels ("+i+","+(j-1)+") and ("+i+","+j+") are of the same color " + newColors[r1] + ".");
                    return -1.0;
                }
            }
        }

        // score needs to store the number of pixels which had to be recolored
        // and the total number of colors used for coloring
        int recolor = 0;
        for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j)
            if (getColor(true, i, j) != getColor(false, i, j))
                recolor++;
        HashSet<Integer> distinct = new HashSet<>();
        for (int i = 0; i < R; ++i)
            distinct.add(newColors[i]);
        if (debug) {
            addFatalError("Number of colors used: " + distinct.size());
            addFatalError("Number of pixels recolored: " + recolor + " / " + (H * W));
        }
        // 100000 * total # of colors + number of pixels recolored
        return 100000.0 * distinct.size() + recolor;
      }
      catch (Exception e) { 
        addFatalError("An exception occurred while trying to get your program's results.");
        e.printStackTrace(); 
        return -1.0;
      }
    }
// ------------- visualization part ------------
    static String exec, fileName;
    static boolean vis, debug;
    static Process proc;
    InputStream is;
    OutputStream os;
    BufferedReader br;
    static int SZ, SZX, SZY;
    // -----------------------------------------
    int[] recolor(int H, int[] regions, int[] oldColors) throws IOException {
        StringBuffer sb = new StringBuffer();
        sb.append(H).append("\n");
        sb.append(regions.length).append("\n");
        for (int i = 0; i < regions.length; ++i) {
            sb.append(regions[i]).append("\n");
        }
        sb.append(oldColors.length).append("\n");
        for (int i = 0; i < oldColors.length; ++i) {
            sb.append(oldColors[i]).append("\n");
        }
        os.write(sb.toString().getBytes());
        os.flush();

        // and get the return value
        int N = Integer.parseInt(br.readLine());
        int[] ret = new int[N];
        for (int i = 0; i < N; i++)
            ret[i] = Integer.parseInt(br.readLine());
        return ret;
    }
    // -----------------------------------------
    java.util.List<Color> palette;
    public void draw(boolean newC) {
        BufferedImage bi = new BufferedImage(SZX, SZY,BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = (Graphics2D)bi.getGraphics();

        // background
        g2.setColor(Color.WHITE);
        g2.fillRect(0, 0, SZX, SZY);

        // colored pixels of the map 
        for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            g2.setColor(palette.get(getColor(newC, i, j)));
            g2.fillRect(j * SZ + 1, i * SZ + 1, SZ - 1, SZ - 1);
        }

        // region indices in each pixel
        g2.setColor(Color.BLACK);
        g2.setFont(new Font("Arial",Font.PLAIN,9));
        FontMetrics fm = g2.getFontMetrics();
        char[] ch;
        for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            ch = ("" + regions[i][j] + (newC ? " / " + getColor(newC, i, j) : "")).toCharArray();
            int h = i * SZ + SZ/2 + fm.getHeight()/2 - 2;
            g2.drawChars(ch, 0, ch.length, j * SZ + 2, h);
        }

        // borders between pixels
        for (int i = 0; i <= H; i++)
            g2.drawLine(0,i*SZ,W*SZ,i*SZ);
        for (int i = 0; i <= W; i++)
            g2.drawLine(i*SZ,0,i*SZ,H*SZ);

        // if two adjacent pixels have a bad border (different regions same color),
        // paint the border bold black, otherwise paint it less bold black
        g2.setStroke(new BasicStroke(3.0f));
        for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            // horizontal lines
            if (i > 0 && regions[i][j] != regions[i-1][j]) {
                int c0 = getColor(newC, i, j);
                int c1 = getColor(newC, i-1, j);
                if (c0 == c1)
                    g2.drawLine(j*SZ, i*SZ, (j+1)*SZ, i*SZ);
            }
            // vertical lines
            if (j > 0 && regions[i][j] != regions[i][j-1]) {
                int c0 = getColor(newC, i, j);
                int c1 = getColor(newC, i, j-1);
                if (c0 == c1)
                    g2.drawLine(j*SZ, i*SZ, j*SZ, (i+1)*SZ);
            }
        }

        try {
            new File("vis").mkdirs();
            ImageIO.write(bi,"png", new File("vis/" + fileName + (newC ? "-res" : "-in") + ".png"));
        } catch (Exception e) { e.printStackTrace(); }
    }
    // -----------------------------------------
    public MapRecoloringVis(String seed) {
      try {
        if (exec != null) {
            try {
                Runtime rt = Runtime.getRuntime();
                proc = rt.exec(exec);
                os = proc.getOutputStream();
                is = proc.getInputStream();
                br = new BufferedReader(new InputStreamReader(is));
                new ErrorReader(proc.getErrorStream()).start();
            } catch (Exception e) { e.printStackTrace(); }
        }
        System.out.println("Score = " + String.format("%.0f", runTest(seed)));
        if (proc != null)
            try { proc.destroy(); } 
            catch (Exception e) { e.printStackTrace(); }
      }
      catch (Exception e) { e.printStackTrace(); }
    }
    // -----------------------------------------
    public static void main(String[] args) {
        String seed = "1";
        vis = true;
        SZ = 20;
        if (seed.equals("1"))
            SZ = 30;
        for (int i = 0; i<args.length; i++)
        {   if (args[i].equals("-seed"))
                seed = args[++i];
            if (args[i].equals("-exec"))
                exec = args[++i];
            if (args[i].equals("-novis"))
                vis = false;
            if (args[i].equals("-size"))
                SZ = Integer.parseInt(args[++i]);
            if (args[i].equals("-debug"))
                debug = true;
        }
        if (exec == null)
            vis = true;
        if (vis)
            fileName = seed;
        MapRecoloringVis f = new MapRecoloringVis(seed);
    }
    // -----------------------------------------
    void addFatalError(String message) {
        System.out.println(message);
    }
}

class ErrorReader extends Thread{
    InputStream error;
    public ErrorReader(InputStream is) {
        error = is;
    }
    public void run() {
        try {
            byte[] ch = new byte[50000];
            int read;
            while ((read = error.read(ch)) > 0)
            {   String s = new String(ch,0,read);
                System.out.print(s);
                System.out.flush();
            }
        } catch(Exception e) { }
    }
}
