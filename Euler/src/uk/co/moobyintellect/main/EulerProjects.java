package uk.co.moobyintellect.main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.InputMismatchException;
import java.util.List;
import java.util.Scanner;
import java.util.stream.Collectors;

public class EulerProjects {

	private static final int[] TEST_BASE = { 2, 7, 61 };
	private static final int NUM_BASES = 3;

	private static String[] ONES = { "zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine" }; // Requires
																														// 0
																														// <=
																														// n
																														// <=
																														// 9
	private static String[] TEENS = { "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen",
			"seventeen", "eighteen", "nineteen" };
	private static String[] TENS = { "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety" };

	public static void findMultiples() {

		try (Scanner scanner = new Scanner(System.in)) {
			System.out.print("Enter number: ");
			List<Integer> numsList = new ArrayList<>();
			int num = scanner.nextInt();
			for (int i = 1; i < num; i++) {
				if ((i % 3 == 0) || (i % 5 == 0)) {
					numsList.add(i);
				}
			}
			System.out.println(numsList.stream().map(i -> i.toString()).collect(Collectors.joining(", ", "{", "}")));
			System.out.println(
					String.format("Total sum: %,d", numsList.stream().mapToInt((x) -> x).summaryStatistics().getSum()));
		} catch (InputMismatchException e ) {
			// TODO: handle exception
		}
	}

	public static void fiboSum() {
		List<Integer> numsList = new ArrayList<>();
		Scanner scanner = new Scanner(System.in);
		System.out.print("Enter highest Fibonacci term: ");
		int highestTerm = scanner.nextInt();
		int currFb = 0;
		numsList.add(0);
		numsList.add(1);
		while (currFb < highestTerm) {
			currFb = numsList.get((numsList.size() - 1)) + (numsList.get((numsList.size() - 2)));
			numsList.add(currFb);
		}
		System.out.println(numsList.stream().map(i -> i.toString()).collect(Collectors.joining(",", "{", "}")));
		System.out.println(String.format("Total Sum: %,d",
				numsList.stream().filter(n -> isEven(n)).mapToInt((x) -> x).summaryStatistics().getSum()));
		scanner.close();
	}

	private static boolean isEven(int num) {
		return ((num % 2) == 0);
	}

	public static void primeFactors() {
		List<Integer> factors = new ArrayList<>();
		Scanner scanner = new Scanner(System.in);
		System.out.print("Enter number to find the prime factors: ");
		long input = scanner.nextLong();
		for (int i = 2; i <= input / i; i++) {
			while (input % i == 0) {
				factors.add(i);
				input /= i;
			}
		}
		if (input > 1) {
			factors.add((int) input);
		}
		System.out.println(factors.stream().map(i -> i.toString()).collect(Collectors.joining(",", "{", "}")));
		scanner.close();
	}

	private static boolean isPalindrome(String s) {
		// if odd length, remove middle character
		if (s.length() % 2 != 0) {
			s = s.substring(0, s.length() / 2) + s.substring((s.length() / 2) + 1);
		}
		return new StringBuilder(s.substring(0, s.length() / 2)).reverse().toString()
				.equals(s.substring(s.length() / 2));
	}

	public static void findPalinromeForThreeDigits() {
		long biggest = 0;
		int a = 999, b = 999;
		int num1 = 0, num2 = 0;
		for (int bb = b; bb > 0; bb--) {
			for (int aa = a; aa > 0; aa--) {
				if (isPalindrome(new Long(aa * bb).toString())) {
					if (aa * bb > biggest) {
						biggest = aa * bb;
						num1 = aa;
						num2 = bb;
					}
				}
			}
		}
		System.out.println(String.format("3-digit numbers: %,d : %,d Palindrome Number: %,d", num1, num2, biggest));
	}

	public static void findNthPrime() {
		System.out.print("Enter Prime term to find: ");
		Scanner scanner = new Scanner(System.in);
		int target = scanner.nextInt();
		long start = System.currentTimeMillis();
		int prime = nthPrime(target);
		long stop = System.currentTimeMillis();
		System.out.println(String.format("Prime No. %8d: %10d\nTime: %8dms", target, prime, (stop - start)));
		scanner.close();
	}

	private static long gcd(long a, long b) {
		while (0 != b) {
			long temp = a;
			a = b;
			b = temp % b;
		}
		return a;
	}

	public static void leastCommonMultiple() {
		System.out.print("Enter  range for which you want to find the Leats Common Multiple: ");
		Scanner scanner = new Scanner(System.in);
		long n = scanner.nextLong();
		long multiple = 1;
		for (long i = 2; i < n; i++) {
			multiple *= i / gcd(i, multiple);
		}
		System.out.println(String.format("Least Common Multiple is: %,d", multiple));
		scanner.close();
	}

	public static void diffOfSquares() {
		System.out.print("Enter number range: ");
		Scanner scanner = new Scanner(System.in);
		long num = scanner.nextLong();
		long sofSq = sumOfSquares(num);
		long sqOfSums = squareofSums(num);
		System.out.println(String.format("Sum of Squares: %10d\nSquare of Sums: %10d\nDifference: %10d", sofSq,
				sqOfSums, (sqOfSums - sofSq)));
		scanner.close();
	}

	private static long sumOfSquares(long a) {
		return ((a * (a + 1) * ((2 * a) + 1)) / 6);
	}

	private static long squareofSums(long a) {
		return (long) Math.pow((a * (a + 1)) / 2, 2);
	}

	private static boolean isPrime(int n) {
		if ((n & 1) == 0)
			return n == 2;
		if (n % 3 == 0)
			return n == 3;
		int step = 4, m;
		boolean onlyTD;
		if (n < 40000) {
			m = (int) Math.sqrt(n) + 1;
			onlyTD = true;
		} else {
			m = 100;
			onlyTD = false;
		}
		for (int i = 5; i < m; step = 6 - step, i += step) {
			if (n % i == 0) {
				return false;
			}
		}
		if (onlyTD) {
			return true;
		}
		long md = n, n1 = n - 1, exp = n - 1;
		int s = 0;
		do {
			++s;
			exp >>= 1;
		} while ((exp & 1) == 0);
		// now n-1 = 2^s * exp
		for (int i = 0; i < NUM_BASES; ++i) {
			long r = modPow(TEST_BASE[i], exp, md);
			if (r == 1)
				continue;
			int j;
			for (j = s; j > 0; --j) {
				if (r == n1)
					break;
				r = (r * r) % md;
			}
			if (j == 0)
				return false;
		}
		return true;
	}

	// calculate base^exponent mod modulus
	// |modulus| must be < 2^31.5 and nonzero, exponent nonnegative
	private static long modPow(long base, long exponent, long modulus) {
		if (exponent == 0)
			return 1;
		if (exponent < 0)
			throw new IllegalArgumentException("Can't handle negative exponents");
		if (modulus < 0)
			modulus = -modulus;
		base %= modulus;
		if (base < 0)
			base += modulus;
		long aux = 1;
		while (exponent > 1) {
			if ((exponent & 1) == 1) {
				aux = (aux * base) % modulus;
			}
			base = (base * base) % modulus;
			exponent >>= 1;
		}
		return (aux * base) % modulus;
	}

	public static int nthPrime(int n) {
		if (n < 2)
			return 2;
		if (n == 2)
			return 3;
		int candidate, count, step = 4;
		for (candidate = 5, count = 2; count < n; step = 6 - step, candidate += step) {
			if (isPrime(candidate)) {
				++count;
			}
		}
		return candidate - step;
	}

	public static void findHighestProduct() {
		String x = "7316717653133062491922511967442657474235534919493496983520312774506326239578318016984801869478851843858615607891129494954595017379583319528532088055111254069874715852386305071569329096329522744304355766896648950445244523161731856403098711121722383113622298934233803081353362766142828064444866452387493035890729629049156044077239071381051585930796086670172427121883998797908792274921901699720888093776657273330010533678812202354218097512545405947522435258490771167055601360483958644670632441572215539753697817977846174064955149290862569321978468622482839722413756570560574902614079729686524145351004748216637048440319989000889524345065854122758866688116427171479924442928230863465674813919123162824586178664583591245665294765456828489128831426076900422421902267105562632111110937054421750694165896040807198403850962455444362981230987879927244284909188845801561660979191338754992005240636899125607176060588611646710940507754100225698315520005593572972571636269561882670428252483600823257530420752963450";
		System.out.print("How many digits do you want to find the highest product for?: ");
		Scanner scanner = new Scanner(System.in);
		int numOfDigits = scanner.nextInt();
		long max = 0;
		for (int i = 0; i <= x.length() - numOfDigits; i++) {
			long p = 1;
			int sup = i + numOfDigits;
			for (int j = i; j < sup; j++) {
				p *= Character.getNumericValue(x.charAt(j));
			}
			if (p > max) {
				max = p;
			}
		}
		System.out.println(String.format("Highest %d digit product is %,d", numOfDigits, max));
		scanner.close();
	}

	public static void pythagoreanTriplet() {
		int a, b, c;
		int stop = 1000;
		for (a = 1; a < stop; a++) {
			for (b = a + 1; b < stop; b++) {
				for (c = b + 1; c < stop; c++) {
					if ((a + b + c) == stop) {
						if ((a * a) + (b * b) == (c * c)) {
							System.out.println(
									String.format("a = %,d\nb = %,d\nc = %,d\na*b*c = %,d", a, b, c, (a * b * c)));
						}
					}
				}
			}
		}
	}

	public static void primesSum() {
		List<Integer> primesList = new ArrayList<>();
		System.out.print("Enter highest Prime: ");
		Scanner scanner = new Scanner(System.in);
		int highestPrime = scanner.nextInt();
		primesList.add(2);
		primesList.add(3);
		int step = 4;
		for (int i = 5; i < highestPrime; step = 6 - step, i += step) {
			if (isPrime(i)) {
				primesList.add(i);
			}
		}
		System.out.println(String.format("Sum of all primes below %,d is: %,d", highestPrime,
				primesList.stream().mapToInt((x) -> x).summaryStatistics().getSum()));
		scanner.close();
	}

	public static void greatest20x20Product() {
		int max = 1;
		final int numberInProd = 4;
		int grid[][] = { { 8, 02, 22, 97, 38, 15, 00, 40, 00, 75, 04, 05, 07, 78, 52, 12, 50, 77, 91, 8 },
				{ 49, 49, 99, 40, 17, 81, 18, 57, 60, 87, 17, 40, 98, 43, 69, 48, 04, 56, 62, 00 },
				{ 81, 49, 31, 73, 55, 79, 14, 29, 93, 71, 40, 67, 53, 88, 30, 03, 49, 13, 36, 65 },
				{ 52, 70, 95, 23, 04, 60, 11, 42, 69, 24, 68, 56, 01, 32, 56, 71, 37, 02, 36, 91 },
				{ 22, 31, 16, 71, 51, 67, 63, 89, 41, 92, 36, 54, 22, 40, 40, 28, 66, 33, 13, 80 },
				{ 24, 47, 32, 60, 99, 03, 45, 02, 44, 75, 33, 53, 78, 36, 84, 20, 35, 17, 12, 50 },
				{ 32, 98, 81, 28, 64, 23, 67, 10, 26, 38, 40, 67, 59, 54, 70, 66, 18, 38, 64, 70 },
				{ 67, 26, 20, 68, 02, 62, 12, 20, 95, 63, 94, 39, 63, 8, 40, 91, 66, 49, 94, 21 },
				{ 24, 55, 58, 05, 66, 73, 99, 26, 97, 17, 78, 78, 96, 83, 14, 88, 34, 89, 63, 72 },
				{ 21, 36, 23, 9, 75, 00, 76, 44, 20, 45, 35, 14, 00, 61, 33, 97, 34, 31, 33, 95 },
				{ 78, 17, 53, 28, 22, 75, 31, 67, 15, 94, 03, 80, 04, 62, 16, 14, 9, 53, 56, 92 },
				{ 16, 39, 05, 42, 96, 35, 31, 47, 55, 58, 88, 24, 00, 17, 54, 24, 36, 29, 85, 57 },
				{ 86, 56, 00, 48, 35, 71, 89, 07, 05, 44, 44, 37, 44, 60, 21, 58, 51, 54, 17, 58 },
				{ 19, 80, 81, 68, 05, 94, 47, 69, 28, 73, 92, 13, 86, 52, 17, 77, 04, 89, 55, 40 },
				{ 04, 52, 8, 83, 97, 35, 99, 16, 07, 97, 57, 32, 16, 26, 26, 79, 33, 27, 98, 66 },
				{ 88, 36, 68, 87, 57, 62, 20, 72, 03, 46, 33, 67, 46, 55, 12, 32, 63, 93, 53, 69 },
				{ 04, 42, 16, 73, 38, 25, 39, 11, 24, 94, 72, 18, 8, 46, 29, 32, 40, 62, 76, 36 },
				{ 20, 69, 36, 41, 72, 30, 23, 88, 34, 62, 99, 69, 82, 67, 59, 85, 74, 04, 36, 16 },
				{ 20, 73, 35, 29, 78, 31, 90, 01, 74, 31, 49, 71, 48, 86, 81, 16, 23, 57, 05, 54 },
				{ 01, 70, 54, 71, 83, 51, 54, 69, 16, 92, 33, 48, 61, 43, 52, 01, 89, 19, 67, 48 } };

		for (int col = 0; col < grid[0].length; col++) {
			for (int row = 0; row < grid[1].length; row++) {

				// Check vertically
				if (row < grid[0].length - numberInProd) {
					int tmp = 1;
					for (int i = 0; i < numberInProd; i++) {
						tmp *= grid[col][row + i];
					}
					max = Math.max(max, tmp);
				}

				// Check horizontally
				if (col < grid[1].length - numberInProd) {
					int tmp = 1;
					for (int i = 0; i < numberInProd; i++) {
						tmp *= grid[col + i][row];
					}
					max = Math.max(max, tmp);
				}

				// Check diagonally upwards / forwards
				if ((col < grid[0].length - numberInProd) && (row >= numberInProd)) {
					int tmp = 1;
					for (int i = 0; i < numberInProd; i++) {
						tmp *= grid[col + i][row - i];
					}
					max = Math.max(max, tmp);
				}

				// Check diagonally downwards / forwards
				if ((row < grid[0].length - numberInProd) && (col < grid[1].length - numberInProd)) {
					int tmp = 1;
					for (int i = 0; i < numberInProd; i++) {
						tmp *= grid[col + i][row + i];
					}
					max = Math.max(max, tmp);
				}
			}
		}
		System.out.println(String.format("Greatest product of four adjacent numbers in grid is %,d", max));
	}

	public static void getDivisorsForTriangular() {
		long begin = System.currentTimeMillis();
		int iter = 1;
		int size = 1;
		int numFactors = 0;
		while (numFactors <= 500) {
			numFactors = countFactors(size);
			iter++;
			size += iter;
		}
		long end = System.currentTimeMillis();
		System.out.println(String.format("First triangle number with over 500 divisors is: %,d\nTime taken: %,dms",
				(size - iter), (end - begin)));
	}

	private static int countFactors(int f) {
		int factors = 0;
		for (int i = 1; i <= Math.sqrt(f); i++) {
			if (f % i == 0) {
				factors += 2;
			}
		}
		return factors;
	}

	public static void first10Digits() {
		List<BigInteger> bigArray = new ArrayList<BigInteger>(
				Arrays.asList(new BigInteger("37107287533902102798797998220837590246510135740250"),
						new BigInteger("46376937677490009712648124896970078050417018260538"),
						new BigInteger("74324986199524741059474233309513058123726617309629"),
						new BigInteger("91942213363574161572522430563301811072406154908250"),
						new BigInteger("23067588207539346171171980310421047513778063246676"),
						new BigInteger("89261670696623633820136378418383684178734361726757"),
						new BigInteger("28112879812849979408065481931592621691275889832738"),
						new BigInteger("44274228917432520321923589422876796487670272189318"),
						new BigInteger("47451445736001306439091167216856844588711603153276"),
						new BigInteger("70386486105843025439939619828917593665686757934951"),
						new BigInteger("62176457141856560629502157223196586755079324193331"),
						new BigInteger("64906352462741904929101432445813822663347944758178"),
						new BigInteger("92575867718337217661963751590579239728245598838407"),
						new BigInteger("58203565325359399008402633568948830189458628227828"),
						new BigInteger("80181199384826282014278194139940567587151170094390"),
						new BigInteger("35398664372827112653829987240784473053190104293586"),
						new BigInteger("86515506006295864861532075273371959191420517255829"),
						new BigInteger("71693888707715466499115593487603532921714970056938"),
						new BigInteger("54370070576826684624621495650076471787294438377604"),
						new BigInteger("53282654108756828443191190634694037855217779295145"),
						new BigInteger("36123272525000296071075082563815656710885258350721"),
						new BigInteger("45876576172410976447339110607218265236877223636045"),
						new BigInteger("17423706905851860660448207621209813287860733969412"),
						new BigInteger("81142660418086830619328460811191061556940512689692"),
						new BigInteger("51934325451728388641918047049293215058642563049483"),
						new BigInteger("62467221648435076201727918039944693004732956340691"),
						new BigInteger("15732444386908125794514089057706229429197107928209"),
						new BigInteger("55037687525678773091862540744969844508330393682126"),
						new BigInteger("18336384825330154686196124348767681297534375946515"),
						new BigInteger("80386287592878490201521685554828717201219257766954"),
						new BigInteger("78182833757993103614740356856449095527097864797581"),
						new BigInteger("16726320100436897842553539920931837441497806860984"),
						new BigInteger("48403098129077791799088218795327364475675590848030"),
						new BigInteger("87086987551392711854517078544161852424320693150332"),
						new BigInteger("59959406895756536782107074926966537676326235447210"),
						new BigInteger("69793950679652694742597709739166693763042633987085"),
						new BigInteger("41052684708299085211399427365734116182760315001271"),
						new BigInteger("65378607361501080857009149939512557028198746004375"),
						new BigInteger("35829035317434717326932123578154982629742552737307"),
						new BigInteger("94953759765105305946966067683156574377167401875275"),
						new BigInteger("88902802571733229619176668713819931811048770190271"),
						new BigInteger("25267680276078003013678680992525463401061632866526"),
						new BigInteger("36270218540497705585629946580636237993140746255962"),
						new BigInteger("24074486908231174977792365466257246923322810917141"),
						new BigInteger("91430288197103288597806669760892938638285025333403"),
						new BigInteger("34413065578016127815921815005561868836468420090470"),
						new BigInteger("23053081172816430487623791969842487255036638784583"),
						new BigInteger("11487696932154902810424020138335124462181441773470"),
						new BigInteger("63783299490636259666498587618221225225512486764533"),
						new BigInteger("67720186971698544312419572409913959008952310058822"),
						new BigInteger("95548255300263520781532296796249481641953868218774"),
						new BigInteger("76085327132285723110424803456124867697064507995236"),
						new BigInteger("37774242535411291684276865538926205024910326572967"),
						new BigInteger("23701913275725675285653248258265463092207058596522"),
						new BigInteger("29798860272258331913126375147341994889534765745501"),
						new BigInteger("18495701454879288984856827726077713721403798879715"),
						new BigInteger("38298203783031473527721580348144513491373226651381"),
						new BigInteger("34829543829199918180278916522431027392251122869539"),
						new BigInteger("40957953066405232632538044100059654939159879593635"),
						new BigInteger("29746152185502371307642255121183693803580388584903"),
						new BigInteger("41698116222072977186158236678424689157993532961922"),
						new BigInteger("62467957194401269043877107275048102390895523597457"),
						new BigInteger("23189706772547915061505504953922979530901129967519"),
						new BigInteger("86188088225875314529584099251203829009407770775672"),
						new BigInteger("11306739708304724483816533873502340845647058077308"),
						new BigInteger("82959174767140363198008187129011875491310547126581"),
						new BigInteger("97623331044818386269515456334926366572897563400500"),
						new BigInteger("42846280183517070527831839425882145521227251250327"),
						new BigInteger("55121603546981200581762165212827652751691296897789"),
						new BigInteger("32238195734329339946437501907836945765883352399886"),
						new BigInteger("75506164965184775180738168837861091527357929701337"),
						new BigInteger("62177842752192623401942399639168044983993173312731"),
						new BigInteger("32924185707147349566916674687634660915035914677504"),
						new BigInteger("99518671430235219628894890102423325116913619626622"),
						new BigInteger("73267460800591547471830798392868535206946944540724"),
						new BigInteger("76841822524674417161514036427982273348055556214818"),
						new BigInteger("97142617910342598647204516893989422179826088076852"),
						new BigInteger("87783646182799346313767754307809363333018982642090"),
						new BigInteger("10848802521674670883215120185883543223812876952786"),
						new BigInteger("71329612474782464538636993009049310363619763878039"),
						new BigInteger("62184073572399794223406235393808339651327408011116"),
						new BigInteger("66627891981488087797941876876144230030984490851411"),
						new BigInteger("60661826293682836764744779239180335110989069790714"),
						new BigInteger("85786944089552990653640447425576083659976645795096"),
						new BigInteger("66024396409905389607120198219976047599490197230297"),
						new BigInteger("64913982680032973156037120041377903785566085089252"),
						new BigInteger("16730939319872750275468906903707539413042652315011"),
						new BigInteger("94809377245048795150954100921645863754710598436791"),
						new BigInteger("78639167021187492431995700641917969777599028300699"),
						new BigInteger("15368713711936614952811305876380278410754449733078"),
						new BigInteger("40789923115535562561142322423255033685442488917353"),
						new BigInteger("44889911501440648020369068063960672322193204149535"),
						new BigInteger("41503128880339536053299340368006977710650566631954"),
						new BigInteger("81234880673210146739058568557934581403627822703280"),
						new BigInteger("82616570773948327592232845941706525094512325230608"),
						new BigInteger("22918802058777319719839450180888072429661980811197"),
						new BigInteger("77158542502016545090413245809786882778948721859617"),
						new BigInteger("72107838435069186155435662884062257473692284509516"),
						new BigInteger("20849603980134001723930671666823555245252804609722"),
						new BigInteger("53503534226472524250874054075591789781264330331690")));
		long begin = System.currentTimeMillis();
		System.out.println(String.format("1st 10 digits of sum: %s",
				bigArray.stream().reduce(BigInteger.ZERO, BigInteger::add).toString().substring(0, 10)));
		long end = System.currentTimeMillis();
		System.out.println(String.format("Time taken: %,dms", end - begin));
	}

	public static void collatzSequence() {
		long begin = System.currentTimeMillis();
		int answer = -1;
		int chain = 0;
		for (int i = 999999; i > 0; i--) {
			int tmpChain = 1;
			long tmp = i;
			while (tmp > 1) {
				if (tmp % 2 == 0) {
					tmp /= 2;
				} else {
					tmp = (tmp * 3) + 1;
				}
				tmpChain++;
			}
			if (tmpChain > chain) {
				chain = tmpChain;
				answer = i;
			}
		}
		long end = System.currentTimeMillis();
		System.out.println(String.format("Collatz Term with longest chain is %,d with %,d terms.\nTime taken: %,dms",
				answer, chain, (end - begin)));
	}

	private static long binomialCoefficient(int n, int k) {
		/*
		 * N-choose-k combinatorics: (n! / (k! * (n - k)!)) where: n is the
		 * number of moves, k is the number of down and right moves required (20
		 * each)
		 */

		if (k > (n - k)) {
			k = n - k;
		}
		long c = 1;
		for (int i = 0; i < k; i++) {
			c = c * (n - i);
			c = c / (i + 1);
		}
		return c;
	}

	public static void findNumberofMoves() {
		long begin = System.currentTimeMillis();
		System.out.println(String.format("Number of moves: %,d", binomialCoefficient(40, 20)));
		long end = System.currentTimeMillis();
		System.out.println(String.format("Time taken: %,dms", (end - begin)));
	}

	private static BigInteger pow(BigInteger base, BigInteger exponent) {
		BigInteger result = BigInteger.ONE;
		while (exponent.signum() > 0) {
			if (exponent.testBit(0)) {
				result = result.multiply(base);
			}
			base = base.multiply(base);
			exponent = exponent.shiftRight(1);
		}
		return result;
	}

	public static void powerDigitSum() {
		long begin = System.currentTimeMillis();
		BigInteger result = BigInteger.ZERO;
		BigInteger two = BigInteger.valueOf(2);
		BigInteger one_thousand = BigInteger.valueOf(1000);

		result = pow(two, one_thousand);

		String string_result = new String(result.toString());

		int answer = 0;
		for (int i = 0; i < string_result.length(); i++) {
			Character c = new Character(string_result.charAt(i));
			String s = c.toString();
			int n = Integer.parseInt(s);
			answer += n;
		}
		long end = System.currentTimeMillis();
		System.out.println(String.format("Sum of digits for 2 to the power 1000 is: %,d\nTime Taken: %,dms", answer,
				(end - begin)));
	}

	public static void numberLetterCounts() {
		long begin = System.currentTimeMillis();
		int sum = 0;
		for (int i = 1; i <= 1000; i++) {
			sum += toEnglish(i).length();
		}
		long end = System.currentTimeMillis();
		System.out.println(
				String.format("Number of Letters used in 1 - 1000 is: %,d\nTime Taken: %,dms", sum, (end - begin)));
	}

	// requires 0 <= n <= 99999
	private static String toEnglish(int n) {
		if (n < 0 || n > 99999) {
			throw new IllegalArgumentException();
		}
		if (n < 100) {
			return tens(n);
		} else {
			String big = "";
			if (n >= 1000) {
				big += tens(n / 1000) + "thousand";
			}
			if (n / 100 % 10 != 0) {
				big += ONES[n / 100 % 10] + "hundred";
			}
			return big + (n % 100 != 0 ? "and" + tens(n % 100) : "");
		}
	}

	// requires 0 <= n ,<= 99
	private static String tens(int n) {
		if (n < 10) {
			return ONES[n];
		} else if (n < 20) { // teens
			return TEENS[n - 10];
		} else {
			return TENS[n / 10 - 2] + (n % 10 != 0 ? ONES[n % 10] : "");
		}
	}

	private static int[][] readTree(int heightOfTree, String fileName) throws IOException {
		int[][] tree = new int[heightOfTree][];
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		for (int i = 0; i < heightOfTree; i++) {
			tree[i] = new int[i + 1];
			String[] values = reader.readLine().split(" ");
			for (int j = 0; j <= i; j++) {
				tree[i][j] = Integer.parseInt(values[j]);
			}
		}
		return tree;
	}

}
