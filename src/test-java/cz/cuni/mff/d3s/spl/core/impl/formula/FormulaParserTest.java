/*
 * Copyright 2013 Charles University in Prague
 * Copyright 2013 Vojtech Horky
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package cz.cuni.mff.d3s.spl.core.impl.formula;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import cz.cuni.mff.d3s.spl.core.Data;
import cz.cuni.mff.d3s.spl.core.Formula;
import cz.cuni.mff.d3s.spl.core.MathematicalInterpretation;
import cz.cuni.mff.d3s.spl.core.Result;
import cz.cuni.mff.d3s.spl.core.impl.DataForTest;
import cz.cuni.mff.d3s.spl.core.impl.InterpretationForTests;
import cz.cuni.mff.d3s.spl.core.impl.formula.SplFormula.SplParseException;

@RunWith(Parameterized.class)
public class FormulaParserTest {
	protected static final int SAMPLE_COUNT_GOOD = InterpretationForTests.MINIMUM_SAMPLES_REQUIRED * 2;
	protected static final int SAMPLE_COUNT_BAD = InterpretationForTests.MINIMUM_SAMPLES_REQUIRED / 2;
	
	/* We use these constants only for nicer formatting. */
	protected static final Result TRU = Result.TRUE;
	protected static final Result FAL = Result.FALSE;
	protected static final Result UNK = Result.CANNOT_COMPUTE;

	@Parameters(name = "{index}: {0} {1}")
	public static Collection<Object[]> createParameters() {
		return Arrays.asList(new Object[][] {
			{ TRU, "low < high" },
			
			{ TRU, "high > low" },
			
			{ TRU, "low < medium && medium < high" },
			
			{ TRU, "(low < medium)" },
			{ TRU, "(high < medium) || (medium < high)" },
			{ TRU, "low < high && (low > medium || high > medium)" },
			
			{ TRU, "low < medium => low < medium" },
			{ FAL, "low < medium => medium < low" },
			{ TRU, "medium < low => low > high" },
		});
	}
	
	private final String formulaAsString;
	private final Result expectedResult;
	
	protected MathematicalInterpretation interpretation;
	protected Data low;
	protected Data medium;
	protected Data high;
	protected Data empty;
	
	public FormulaParserTest(final Result result, final String formula) {
		formulaAsString = formula;
		expectedResult = result;
	}
	
	@Before
	public void setupInterpretation() {
		interpretation = new InterpretationForTests();
	}
	
	@Before
	public void setupData() {
		low = new DataForTest(10, SAMPLE_COUNT_GOOD);
		medium = new DataForTest(20, SAMPLE_COUNT_GOOD);
		high = new DataForTest(30, SAMPLE_COUNT_GOOD);
		empty = new DataForTest(0, SAMPLE_COUNT_BAD);
	}
	
	@Test
	public void parse() throws SplParseException {
		Formula formula = SplFormula.create(formulaAsString);
		
		assertNotNull(formula);
		
		/* First, set the proper interpretation. */
		formula.setInterpreation(new InterpretationForTests());

		/* Bind the provided data sources. */
		bindIfPresent(formula, formulaAsString);
		
		/* And evaluate. */
		assertEquals(expectedResult, formula.evaluate());
	}
	
	private void bindIfPresent(Formula parsed, String formula) {
		bindIfPresent(parsed, formula, low, "low");
		bindIfPresent(parsed, formula, medium, "medium");
		bindIfPresent(parsed, formula, high, "high");
		bindIfPresent(parsed, formula, empty, "empty");
	}
	
	private void bindIfPresent(Formula parsed, String formula, Data source, String variable) {
		if (formula.contains(variable)) {
			parsed.bind(variable, source);
		}
	}
}
