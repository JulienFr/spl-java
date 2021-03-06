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
package cz.cuni.mff.d3s.spl.instrumentation;

public class ExtraArgument {
	private enum Kind {
		NULL,
		THIS,
		FIELD,
		PARAMETER,
	};
	private final Kind kind;
	private final String name;
	private final int index;
	
	public static final ExtraArgument PASS_THIS = new ExtraArgument(Kind.THIS, null, -1);
	public static final ExtraArgument PASS_NULL = new ExtraArgument(Kind.NULL, null, -1);
	public static final ExtraArgument PASS_THROUGH_PARAMETER_1 = new ExtraArgument(Kind.PARAMETER, null, 1);
	public static final ExtraArgument PASS_THROUGH_PARAMETER_2 = new ExtraArgument(Kind.PARAMETER, null, 2);
	public static final ExtraArgument PASS_THROUGH_PARAMETER_3 = new ExtraArgument(Kind.PARAMETER, null, 3);
	public static final ExtraArgument PASS_THROUGH_PARAMETER_4 = new ExtraArgument(Kind.PARAMETER, null, 4);
	
	private ExtraArgument(Kind kind, String name, int index) {
		this.kind = kind;
		this.name = name;
		this.index = index;
	}
	
	public void accept(ExtraArgumentVisitor visitor) {
		switch (kind) {
		case NULL:
			visitor.visitNull();
			return;
		case THIS:
			visitor.visitThis();
			return;
		case FIELD:
			visitor.visitField(name);
			return;
		case PARAMETER:
			visitor.visitParameter(index);
			return;
		}
	}
	
	public static ExtraArgument createField(String name) {
		return new ExtraArgument(Kind.FIELD, name, -1);
	}
	
	public static ExtraArgument createParameter(int position) {
		return new ExtraArgument(Kind.PARAMETER, null, position);
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + index;
		result = prime * result + ((kind == null) ? 0 : kind.hashCode());
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (!(obj instanceof ExtraArgument)) {
			return false;
		}
		ExtraArgument other = (ExtraArgument) obj;
		if (index != other.index) {
			return false;
		}
		if (kind != other.kind) {
			return false;
		}
		if (name == null) {
			if (other.name != null) {
				return false;
			}
		} else if (!name.equals(other.name)) {
			return false;
		}
		return true;
	}
}
