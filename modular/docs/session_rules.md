# Session Rules

Enforce strictly in every session.

1. **Instrument before reasoning.** When debugging, add print/assert statements first, run the program, read the output, *then* reason about the fix. Never spend more than one round-trip reasoning about a bug without observed evidence.

2. **Incremental output.** After every non-trivial tool call (compile, run, edit), stop and report the result in ≤5 lines before continuing. Do not chain more than ~3 tool calls silently.

3. **Checkpoint before large writes.** Before writing more than ~50 lines of new code, state what you are about to do in one sentence and wait. Do not produce multiple large files in a single response.

4. **Fail fast.** If a step produces unexpected output, stop immediately, report it, and ask before proceeding.

5. **No speculative work.** Do not implement the next layer while the current one has unverified behaviour.
