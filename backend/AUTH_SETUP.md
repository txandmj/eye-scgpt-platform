# Request User Binding (Google headers + test fallback)

This layer binds each request to a user row in the `users` table using headers. If headers are missing (local/dev), a single test user is created/reused.

## Headers
- `X-Google-Sub` (string): unique subject identifier from Google. If present, it is used as the unique key.
- `X-User-Email` (optional): user email.
- `X-User-Name` (optional): display name.

## Behavior
- If `X-Google-Sub` is present: upsert into `users` by `google_sub` and return that user.
- If missing: use a single test user persisted once:
  - `google_sub = test-sub`
  - `email = test@example.com`
  - `display_name = Test User`

The dependency returns a lightweight dict: `{ id, google_sub, email, display_name }`.

## Example (curl)

With real headers:
```
curl -H "X-Google-Sub: 1234567890abcdef" \
     -H "X-User-Email: alice@example.com" \
     -H "X-User-Name: Alice" \
     http://localhost:8000/api/jobs
```

Without headers (local dev fallback):
```
curl http://localhost:8000/api/jobs
```

## Logging
Each request logs which user was bound (either Google user or the fallback test user).

## Notes
- No API response shapes are changed; existing tests should continue to pass.
- To switch to PostgreSQL, set `DATABASE_URL` accordingly; otherwise SQLite is used by default.


