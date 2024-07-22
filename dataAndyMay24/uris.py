import uuid

def generate_mini_uri(length=8):
    # Generate a UUID
    hash_code = str(uuid.uuid4()).replace('-', '')
    # Extract a substring of the desired length
    return hash_code[:length]

# Generate 10 mini URIs
mini_uris = [generate_mini_uri() for _ in range(10)]

for i, mini_uri in enumerate(mini_uris, 1):
    print(mini_uri)
