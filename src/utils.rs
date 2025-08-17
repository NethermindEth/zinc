/// Reinterpret a `Vec<Self>` as a `Vec<Target>` **without copying**.
///
/// This trait provides a zero-copy conversion between vectors whose element
/// types are *layout-compatible*. It is meant for cases like transparent
/// newtypes (`#[repr(transparent)] struct Newtype(Inner);`) where the two
/// types have identical ABI (size and alignment) and the *value validity*
/// requirements of `Target` are **no stricter** than those of `Self`.
///
/// # What this does
///
/// Given ownership of a `Vec<Self>`, the implementation:
/// - Prevents the original vector from being dropped (`ManuallyDrop`),
/// - Reuses its allocation, pointer, length, and capacity unchanged,
/// - Reinterprets the element pointer as `*mut Target`,
/// - Constructs a new `Vec<Target>` with the same buffer.
///
/// No bytes are moved or rewritten. This is O(1).
///
/// # Why this trait is `unsafe`
///
/// This is an **`unsafe trait`** because **implementing** it asserts a
/// *per-type* contract that the compiler cannot verify. A bad impl can cause
/// undefined behavior in safe code that *uses* the impl. See
/// the *Implementor requirements* section below.
///
/// Additionally, the **method** is `unsafe` because each **call site** must
/// ensure the conversion is valid for the particular `Vec<Self>` value being
/// reinterpreted. (You may choose to make the method safe and keep only the
/// trait unsafe; keeping both is conservative but redundant.)
///
/// # Safety
///
/// ## Caller requirements (for each call to [`reinterpret_vector`])
///
/// You must uphold all of the following at the call site:
///
/// 1. **Layout equivalence**
///    `size_of::<Self>() == size_of::<Target>()` and
///    `align_of::<Self>() == align_of::<Target>()`.
///    This implementation *asserts* these at compile time, failing the build
///    if they do not hold. (Do not remove those checks.)
///
/// 2. **Validity subset**
///    Every bit-pattern that is a valid `Self` value must also be a valid
///    `Target` value. In other words, `Target` must **not** impose stricter
///    validity or niche invariants than `Self`.
///
///    **Counterexamples (do NOT reinterpret into these):**
///    - `*mut T` → `NonNull<T>` (null is valid for `*mut T` but **invalid** for `NonNull<T>`).
///    - `u64` → `core::num::NonZeroU64`.
///    - `Option<NonNull<T>>` has niche optimization differences compared to raw pointers.
///
/// 3. **Drop semantics remain sound**
///    Dropping `Vec<Target>` must be sound for all elements coming from
///    `Vec<Self>`. If `Target` has a custom `Drop` impl that relies on
///    invariants not guaranteed by `Self`, the conversion is unsound.
///
/// 4. **Ownership/uniqueness is preserved**
///    After calling this function, the original `Vec<Self>` **must not** be
///    used again. (The function takes it by value and places it in
///    `ManuallyDrop` to prevent double-free.)
///
/// ## Implementor requirements (for each `unsafe impl`)
///
/// When you write `unsafe impl ReinterpretVector<Target> for Self`, you are
/// promising that for **any** `Vec<Self>`:
///
/// - `Self` and `Target` are layout-identical (same size and alignment) —
///   e.g. `Target` is `#[repr(transparent)]` over `Self`, or is otherwise
///   guaranteed to share ABI with `Self`.
/// - The set of valid `Self` values is a subset of the valid `Target` values.
///   No stricter validity or niche invariants are introduced by `Target`.
/// - Dropping a `Target` produced from arbitrary valid `Self` bytes is sound
///   (no additional assumptions in `Drop`).
///
/// If you cannot *globally* guarantee these at the type level, do **not**
/// implement this trait. Prefer a dedicated constructor that validates data or
/// a fallible conversion (`TryFrom<Vec<Self>> for Vec<Target>` that maps
/// elements).
///
/// ### Recommended: seal the trait
///
/// To prevent downstream crates from writing unsound impls, consider “sealing”
/// the trait (making it unimplementable outside your crate) and providing only
/// vetted `unsafe impl`s for your own types.
///
/// # Edge cases and notes
///
/// - **Zero-Sized Types (ZSTs):** Reinterpreting between ZSTs is fine if
///   alignment matches. `Vec<T>` for ZST `T` uses a dangling non-null pointer;
///   `Vec::from_raw_parts` supports this. The compile-time size/alignment
///   asserts still gate the conversion.
/// - **Uninhabited types (`!`):** `Vec<!>` cannot exist; do not attempt.
/// - **Provenance / aliasing:** We preserve the original allocation, pointer,
///   length, and capacity. Only the *type* of the element pointer changes
///   (with equal alignment), which is allowed under Rust’s aliasing rules for
///   layout-identical types when validity is preserved.
/// - **Send/Sync auto-traits:** `Vec<Target>` may have different `Send`/`Sync`
///   auto-trait behavior than `Vec<Self>`. That is a *type-system* effect at
///   compile time; the runtime conversion does not change the buffer contents.
/// - **Alternative (borrowed view):** If you only need a temporary view,
///   prefer reinterpreting slices to avoid ownership transfer:
///   ```ignore
///   // SAFETY: same layout + lifetime/aliasing upheld
///   let view: &[Target] = unsafe {
///       core::slice::from_raw_parts(vec_self.as_ptr().cast::<Target>(), vec_self.len())
///   };
///   ```
///
/// # Examples
///
/// Transparent newtype over a primitive:
///
/// # Implementation notes
///
/// - The internal `const` block uses compile-time `assert!`s to fail builds
///   when layout differs. Keep these assertions to defend against accidental
///   misuse.
/// - The method uses `ManuallyDrop` to avoid dropping `source` before we
///   reconstruct a `Vec<Target>` from its raw parts. This prevents double-free.
/// - `Vec::from_raw_parts` is used in an `unsafe` block because we must
///   guarantee that `(ptr, len, cap)` originated from a `Vec<Target>`-compatible
///   allocation. The trait’s contract + layout asserts provide that guarantee.
/// - Consider enabling `#![deny(unsafe_op_in_unsafe_fn)]` at the crate root so
///   that all unsafe operations within this `unsafe fn` remain explicitly
///   marked and audited.
///
/// # When to prefer a safe, copying conversion
///
/// If `Target` applies *any* additional invariants (non-zero, non-null, range
/// restrictions, encoding constraints, etc.), use a per-element conversion
/// (`map(Into::into)` or `try_map`) that can validate or fail gracefully. Zero-copy
/// reinterprets are only appropriate when the types are truly the same at the
/// byte level.
///
/// ---
pub(crate) unsafe trait ReinterpretVector<Target: Sized>: Sized {
    /// Reinterpret a `Vec<Self>` as a `Vec<Target>` **without copying**.
    ///
    /// See the trait-level documentation for full safety requirements and examples.
    ///
    /// # Safety
    ///
    /// Caller must ensure layout equivalence, validity preservation, and sound
    /// drop semantics for `Target` given elements produced from `Self`. After
    /// calling, the `source` vector must not be used again.
    unsafe fn reinterpret_vector(source: Vec<Self>) -> Vec<Target> {
        use std::mem::ManuallyDrop;
        const {
            // Keep these to enforce layout at compile time.
            assert!(std::mem::size_of::<Self>() == std::mem::size_of::<Target>());
            assert!(std::mem::align_of::<Self>() == std::mem::align_of::<Target>());
        }
        // Prevent `source` from being dropped while we steal its parts.
        let mut v = ManuallyDrop::new(source);
        let (ptr, len, cap) = (v.as_mut_ptr(), v.len(), v.capacity());
        // SAFETY:
        // - ptr/len/cap come from a `Vec<Self>` allocation of `cap` elements,
        // - size/align equality ensures `Vec<Target>` layout compatibility,
        // - ownership of the buffer is transferred, `source` will not be used again.
        unsafe { Vec::from_raw_parts(ptr.cast::<Target>(), len, cap) }
    }

    /// Reinterpret a borrowed slice `&[Self]` as `&[Target]` **without copying**.
    ///
    /// This is an O(1) borrowed view that preserves the original lifetime and
    /// does not take ownership of the buffer.
    ///
    /// # Safety
    ///
    /// The caller must ensure:
    /// 1. **Layout equivalence**: `size_of::<Self>() == size_of::<Target>()` and
    ///    `align_of::<Self>() == align_of::<Target>()`. (Enforced via compile-time asserts.)
    /// 2. **Validity subset**: Every valid `Self` bit-pattern is also a valid `Target`.
    ///    (No stricter validity or niche invariants are introduced by `Target`.)
    /// 3. **Aliasing rules**: The returned `&[Target]` must not be used to violate Rust’s
    ///    aliasing model (e.g., you must not concurrently mutate the same memory via `Self`
    ///    references while holding this `&[Target]`).
    ///
    /// Unlike the vector method, there are no drop-semantics concerns here because no
    /// ownership is transferred; this is a view only.
    unsafe fn reinterpret_slice(source: &[Self]) -> &[Target] {
        const {
            assert!(std::mem::size_of::<Self>() == std::mem::size_of::<Target>());
            assert!(std::mem::align_of::<Self>() == std::mem::align_of::<Target>());
        }
        // SAFETY:
        // - `source.as_ptr()` points to `len` contiguous `Self` elements.
        // - Equal size/align means it is valid to view the same bytes as `Target`.
        // - Lifetime `'a` is preserved; no ownership is taken.
        unsafe { core::slice::from_raw_parts(source.as_ptr().cast::<Target>(), source.len()) }
    }
}
