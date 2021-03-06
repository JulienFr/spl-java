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
package cz.cuni.mff.d3s.spl.annotations;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/** Marks method to be executed at program termination.
 * 
 * The annotated method is expected to be of public static void type
 * with no arguments.
 * <p>
 * The method is registered from the agent.
 * Thus the code will typically run in a slightly different environment
 * than the normal code.
 * <p>
 * The agent uses Runtime.addShutdownHook to register the user method.
 * Abnormal termination of the JVM may cause that the user code is not
 * executed at all.
 * 
 * @see cz.cuni.mff.d3s.spl.annotations
 */
@Target(ElementType.METHOD)
@Retention(RetentionPolicy.RUNTIME)
public @interface AtExit {
	
}
